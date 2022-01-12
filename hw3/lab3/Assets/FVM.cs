using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;

public class FVM : MonoBehaviour
{
	float dt 			= 0.003f;
    float mass 			= 1;
	float stiffness_0	= 20000.0f;
    float stiffness_1 	= 5000.0f;
    float damp			= 0.999f;

	int[] 		Tet;
	int tet_number;			//The number of tetrahedra

	Vector3[] 	Force;
	Vector3[] 	V;
	Vector3[] 	X;
	int number;				//The number of vertices

	Matrix4x4[] inv_Dm;
    float[] det_Dm;

    //For Laplacian smoothing.
    Vector3[]   V_sum;
	int[]		V_num;

	SVD svd = new SVD();

    private Transform floor;

    HashSet<int>[] vertices_graph;

    // Start is called before the first frame update
    void Start()
    {
    	// FILO IO: Read the house model from files.
    	// The model is from Jonathan Schewchuk's Stellar lib.
    	{
    		string fileContent = File.ReadAllText("Assets/house2.ele");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		
    		tet_number=int.Parse(Strings[0]);
        	Tet = new int[tet_number*4];

    		for(int tet=0; tet<tet_number; tet++)
    		{
				Tet[tet*4+0]=int.Parse(Strings[tet*5+4])-1;
				Tet[tet*4+1]=int.Parse(Strings[tet*5+5])-1;
				Tet[tet*4+2]=int.Parse(Strings[tet*5+6])-1;
				Tet[tet*4+3]=int.Parse(Strings[tet*5+7])-1;
			}
    	}
    	{
			string fileContent = File.ReadAllText("Assets/house2.node");
    		string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
    		number = int.Parse(Strings[0]);
    		X = new Vector3[number];
       		for(int i=0; i<number; i++)
       		{
       			X[i].x=float.Parse(Strings[i*5+5])*0.4f;
       			X[i].y=float.Parse(Strings[i*5+6])*0.4f;
       			X[i].z=float.Parse(Strings[i*5+7])*0.4f;
       		}
    		//Centralize the model.
	    	Vector3 center=Vector3.zero;
	    	for(int i=0; i<number; i++)		center+=X[i];
	    	center=center/number;
	    	for(int i=0; i<number; i++)
	    	{
	    		X[i]-=center;
	    		float temp=X[i].y;
	    		X[i].y=X[i].z;
	    		X[i].z=temp;
	    	}
		}
        /*tet_number=1;
        Tet = new int[tet_number*4];
        Tet[0]=0;
        Tet[1]=1;
        Tet[2]=2;
        Tet[3]=3;

        number=4;
        X = new Vector3[number];
        V = new Vector3[number];
        Force = new Vector3[number];
        X[0]= new Vector3(0, 0, 0);
        X[1]= new Vector3(1, 0, 0);
        X[2]= new Vector3(0, 1, 0);
        X[3]= new Vector3(0, 0, 1);*/

        //Create triangle mesh.
       	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];

        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];

        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }

        int[] triangles = new int[tet_number*12];
        for(int t=0; t<tet_number*4; t++)
        {
        	triangles[t*3+0]=t*3+0;
        	triangles[t*3+1]=t*3+1;
        	triangles[t*3+2]=t*3+2;
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.triangles = triangles;
		mesh.RecalculateNormals ();


		V 	  = new Vector3[number];
        Force = new Vector3[number];
        V_sum = new Vector3[number];
        V_num = new int[number];
        vertices_graph = new HashSet<int>[number];

        if(GameObject.Find("Floor") == null){
            Debug.LogError("Floor doesn't exist in the scene");
        }
        floor = GameObject.Find("Floor").transform;

        //TODO: Need to allocate and assign inv_Dm
        inv_Dm = new Matrix4x4[tet_number];
        det_Dm = new float[tet_number];
        for (int tet = 0; tet < tet_number;++tet){
            Matrix4x4 Dm = Build_Edge_Matrix(tet);

            det_Dm[tet] = Dm[0, 0] * Dm[1, 1] * Dm[2, 2]
                + Dm[0, 1] * Dm[1, 2] * Dm[2, 0]
                + Dm[1, 0] * Dm[2, 1] * Dm[0, 2]
                - Dm[2, 0] * Dm[1, 1] * Dm[0, 2]
                - Dm[0, 1] * Dm[1, 0] * Dm[2, 2]
                - Dm[0, 0] * Dm[1, 2] * Dm[2, 1];
				
            Matrix4x4 adjointM = Matrix4x4.zero;
            adjointM[0, 0] = Dm[1, 1] * Dm[2, 2] - Dm[2, 1] * Dm[1, 2];
            adjointM[1, 0] = -(Dm[1, 0] * Dm[2, 2] - Dm[2, 0] * Dm[1, 2]);
            adjointM[2, 0] = Dm[1, 0] * Dm[2, 1] - Dm[1, 1] * Dm[2, 0];
            
			adjointM[0, 1] = -(Dm[0, 1] * Dm[2, 2] - Dm[2, 1] * Dm[0, 2]);
            adjointM[1, 1] = Dm[0, 0] * Dm[2, 2] - Dm[2, 0] * Dm[0, 2];
            adjointM[2, 1] = -(Dm[0, 0] * Dm[2, 1] - Dm[2, 0] * Dm[0, 1]);

            adjointM[0, 2] = Dm[0, 1] * Dm[1, 2] - Dm[1, 1] * Dm[0, 2];
            adjointM[1, 2] = -(Dm[0, 0] * Dm[1, 2] - Dm[1, 0] * Dm[0, 2]);
            adjointM[2, 2] = Dm[0, 0] * Dm[1, 1] - Dm[1, 0] * Dm[0, 1];

            inv_Dm[tet] = multi(1.0f / det_Dm[tet], adjointM);
            int v0 = Tet[tet * 4 + 0], v1 = Tet[tet * 4 + 1], v2 = Tet[tet * 4 + 2], v3 = Tet[tet * 4 + 3];
            if(vertices_graph[v0] == null){
                vertices_graph[v0] = new HashSet<int>();
            }
            if(vertices_graph[v1] == null){
                vertices_graph[v1] = new HashSet<int>();
            }
            if (vertices_graph[v2] == null)
            {
                vertices_graph[v2] = new HashSet<int>();
            }
            if (vertices_graph[v3] == null)
            {
                vertices_graph[v3] = new HashSet<int>();
            }
            vertices_graph[v0].Add(v1);
            vertices_graph[v0].Add(v2);
            vertices_graph[v0].Add(v3);

            vertices_graph[v1].Add(v0);
            vertices_graph[v1].Add(v2);
            vertices_graph[v1].Add(v3);

            vertices_graph[v2].Add(v0);
            vertices_graph[v2].Add(v1);
            vertices_graph[v2].Add(v3);

            vertices_graph[v3].Add(v0);
            vertices_graph[v3].Add(v1);
            vertices_graph[v3].Add(v2);
        }
    }

    Matrix4x4 Build_Edge_Matrix(int tet)
    {
    	Matrix4x4 ret = Matrix4x4.zero;
        //TODO: Need to build edge matrix here.
        Vector3 X10 = X[Tet[tet * 4 + 0]] - X[Tet[tet * 4 + 1]];
        Vector3 X20 = X[Tet[tet * 4 + 0]] - X[Tet[tet * 4 + 2]];
        Vector3 X30 = X[Tet[tet * 4 + 0]] - X[Tet[tet * 4 + 3]];
        ret[0, 0] = X10.x; ret[1, 0] = X10.y; ret[2, 0] = X10.z;
        ret[0, 1] = X20.x; ret[1, 1] = X20.y; ret[2, 1] = X20.z;
        ret[0, 2] = X30.x; ret[1, 2] = X30.y; ret[2, 2] = X30.z;
        return ret;
    }

    Matrix4x4 subtract(Matrix4x4 a, Matrix4x4 b)
    {
        Matrix4x4 res = Matrix4x4.zero;
        res[0, 0] = a[0, 0] - b[0, 0];
        res[0, 1] = a[0, 1] - b[0, 1];
        res[0, 2] = a[0, 2] - b[0, 2];

        res[1, 0] = a[1, 0] - b[1, 0];
        res[1, 1] = a[1, 1] - b[1, 1];
        res[1, 2] = a[1, 2] - b[1, 2];

        res[2, 0] = a[2, 0] - b[2, 0];
        res[2, 1] = a[2, 1] - b[2, 1];
        res[2, 2] = a[2, 2] - b[2, 2];

        return res;
    }
    Matrix4x4 add(Matrix4x4 a, Matrix4x4 b)
    {
        Matrix4x4 res = Matrix4x4.zero;
        res[0, 0] = a[0, 0] + b[0, 0];
        res[0, 1] = a[0, 1] + b[0, 1];
        res[0, 2] = a[0, 2] + b[0, 2];

        res[1, 0] = a[1, 0] + b[1, 0];
        res[1, 1] = a[1, 1] + b[1, 1];
        res[1, 2] = a[1, 2] + b[1, 2];

        res[2, 0] = a[2, 0] + b[2, 0];
        res[2, 1] = a[2, 1] + b[2, 1];
        res[2, 2] = a[2, 2] + b[2, 2];

        return res;
    }
	Matrix4x4 multi(float coef, Matrix4x4 m){
        Matrix4x4 ret = Matrix4x4.zero;
        ret[0, 0] = coef * m[0, 0];
        ret[0, 1] = coef * m[0, 1];
        ret[0, 2] = coef * m[0, 2];

        ret[1, 0] = coef * m[1, 0];
        ret[1, 1] = coef * m[1, 1];
        ret[1, 2] = coef * m[1, 2];

        ret[2, 0] = coef * m[2, 0];
        ret[2, 1] = coef * m[2, 1];
        ret[2, 2] = coef * m[2, 2];

        return ret;
    }

	Matrix4x4 matrixDot(Matrix4x4 a, Matrix4x4 b){
        Matrix4x4 ret = Matrix4x4.zero;
        for (int i = 0; i < 3;++i){
            for (int j = 0; j < 3;++j){
                for (int k = 0; k < 3;++k){
                    ret[i, j] += a[i, k] * b[k, j];
                }
            }

        }
        return ret;
    }
	
	float trace(Matrix4x4 m){
        return m[0, 0] + m[1, 1] + m[2, 2];
    }

    float massInv = 1.0f;
    Vector3 floorN = Vector3.up;
    float restitution = 0.5f;                 // for collision

    public bool StVK = false;
    void _Update()
    {
    	// Jump up.
		if(Input.GetKeyDown(KeyCode.Space))
    	{
    		for(int i=0; i<number; i++)
    			V[i].y+=0.2f;
    	}

    	for(int i=0; i < number; i++)
    	{
            //TODO: Add gravity to Force.
            Force[i].x = Force[i].z = 0.0f;
            Force[i].y = -9.8f;

            V_sum[i].x = V_sum[i].y = V_sum[i].z = 0;
            V_num[i] = 0;
        }
        float inv_6 = 1.0f / 6.0f;
        if (!StVK)
        {
            for(int tet=0; tet < tet_number; ++tet)
            {
                //TODO: Deformation Gradient
                Matrix4x4 part1 = Build_Edge_Matrix(tet);
                Matrix4x4 F = part1 * inv_Dm[tet];

                //TODO: Green Strain
                Matrix4x4 G = subtract(F.transpose * F, Matrix4x4.identity);

                //TODO: First PK Stress
                Matrix4x4 S = add(multi(stiffness_1, G), multi(stiffness_0 * trace(G) * 0.5f, Matrix4x4.identity));
                Matrix4x4 P = F * S;

                //TODO: Elastic Force
                Matrix4x4 EF = multi(-det_Dm[tet] * inv_6, P * inv_Dm[tet].transpose);
                Vector3 f1 = new Vector3(EF[0, 0], EF[1, 0], EF[2, 0]);
                Vector3 f2 = new Vector3(EF[0, 1], EF[1, 1], EF[2, 1]);
                Vector3 f3 = new Vector3(EF[0, 2], EF[1, 2], EF[2, 2]);
                int v0 = Tet[tet * 4 + 0], v1 = Tet[tet * 4 + 1], v2 = Tet[tet * 4 + 2], v3 = Tet[tet * 4 + 3];
                Force[v1] += f1;
                Force[v2] += f2;
                Force[v3] += f3;
                Force[v0] += -(f1 + f2 + f3);
            }
        }
        else
        {
            for (int tet = 0; tet < tet_number; ++tet)
            {
                // Deformation Gradient
                Matrix4x4 part1 = Build_Edge_Matrix(tet);
                Matrix4x4 F = part1 * inv_Dm[tet];

                // Priciple streches
                Matrix4x4 U = Matrix4x4.zero, S = Matrix4x4.zero, VT = Matrix4x4.zero;
                svd.svd(F, ref U, ref S, ref VT);

                float lambda0 = S[0, 0], lambda1 = S[1, 1], lambda2 = S[2, 2];
                float sqrLambda0 = lambda0 * lambda0
                    , sqrLambda1 = lambda1 * lambda1
                    , sqrLambda2 = lambda2 * lambda2;
                    
                float IC = sqrLambda0 + sqrLambda1 + sqrLambda2;
                // float IIC = sqrLambda0 * sqrLambda1 + sqrLambda0 * sqrLambda2 + sqrLambda1 * sqrLambda2;
                // First PK Stress by StVK model
                // float W = 0.5f * stiffness_0 * (IC - 3) + 0.25f * stiffness_1 * (IIC - 2 * IC + 3);
                float dWd0 = 0.5f * stiffness_0 * lambda0 * (IC - 3) + 0.5f * stiffness_1 * lambda0 * (sqrLambda1 + sqrLambda2 - 2);
                float dWd1 = 0.5f * stiffness_0 * lambda1 * (IC - 3) + 0.5f * stiffness_1 * lambda1 * (sqrLambda0 + sqrLambda2 - 2);
                float dWd2 = 0.5f * stiffness_0 * lambda2 * (IC - 3) + 0.5f * stiffness_1 * lambda2 * (sqrLambda0 + sqrLambda1 - 2);

                // float dWd0 = 0.5f * stiffness_0 * lambda0 * (IC - 3) + stiffness_1 * lambda0 * (sqrLambda0 - 1);
                // float dWd1 = 0.5f * stiffness_0 * lambda1 * (IC - 3) + stiffness_1 * lambda1 * (sqrLambda1 - 1);
                // float dWd2 = 0.5f * stiffness_0 * lambda2 * (IC - 3) + stiffness_1 * lambda2 * (sqrLambda2 - 1);
                Matrix4x4 M = Matrix4x4.zero;
                M[0, 0] = dWd0;
                M[1, 1] = dWd1;
                M[2, 2] = dWd2;
                Matrix4x4 P = U * M * VT;

                // Elastic Force
                Matrix4x4 EF = multi(-det_Dm[tet] * inv_6, P * inv_Dm[tet].transpose);
                Vector3 f1 = new Vector3(EF[0, 0], EF[1, 0], EF[2, 0]);
                Vector3 f2 = new Vector3(EF[0, 1], EF[1, 1], EF[2, 1]);
                Vector3 f3 = new Vector3(EF[0, 2], EF[1, 2], EF[2, 2]);
                int v0 = Tet[tet * 4 + 0], v1 = Tet[tet * 4 + 1], v2 = Tet[tet * 4 + 2], v3 = Tet[tet * 4 + 3];
                Force[v1] += f1;
                Force[v2] += f2;
                Force[v3] += f3;
                Force[v0] += -(f1 + f2 + f3);
            }
        }

        // Laplacian smoothing
        for(int i = 0;i<number;++i){
            foreach(var val in vertices_graph[i]){
                V_sum[i] += V[val];
            }
            V[i] = (V_sum[i] / vertices_graph[i].Count) * 0.2f + V[i] * 0.8f;
        }

        for(int i=0; i<number; i++)
    	{
            //TODO: Update X and V here.
            V[i] += Force[i] * massInv * dt;
            V[i] *= damp;
            X[i] += V[i] * dt;
        
            //TODO: (Particle) collision with floor.
            if(X[i].y < floor.position.y){
                X[i].y = floor.position.y + 0.01f;
                // compute
                Vector3 vn = Vector3.Dot(V[i], floorN) * floorN;
                Vector3 vt = V[i] - vn;
                float a = Mathf.Max(0, 1 - restitution * (1 + restitution) * vn.magnitude / vt.magnitude);
                vn = -restitution * vn;
                vt = a * vt;
                Vector3 vi_new = vn + vt;
                V[i] = vi_new;
            }
        }
    }

    // Update is called once per frame
    void Update()
    {
    	for(int l=0; l<5; l++)
    		 _Update();

    	// Dump the vertex array for rendering.
    	Vector3[] vertices = new Vector3[tet_number*12];
        int vertex_number=0;
        for(int tet=0; tet<tet_number; tet++)
        {
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+0]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        	vertices[vertex_number++]=X[Tet[tet*4+1]];
        	vertices[vertex_number++]=X[Tet[tet*4+2]];
        	vertices[vertex_number++]=X[Tet[tet*4+3]];
        }
        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.vertices  = vertices;
		mesh.RecalculateNormals();
    }
}
