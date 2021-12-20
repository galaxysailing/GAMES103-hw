using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class implicit_model : MonoBehaviour
{
	float 		t 		= 0.0333f;
	float 		mass	= 1;
	float		damping	= 0.99f;
	float 		rho		= 0.8f;
	float 		spring_k = 8000;
	int[] 		E;
	float[] 	L;
	Vector3[] 	V;

    GameObject sphere;

    int n = 21;
    // Start is called before the first frame update
    void Start()
    {
        sphere = GameObject.Find("Sphere");
		if(sphere == null){
            Debug.LogError("Sphere is not in the scene.");
        }

        Mesh mesh = GetComponent<MeshFilter> ().mesh;
		//Resize the mesh.
		Vector3[] X  	= new Vector3[n*n];
		Vector2[] UV 	= new Vector2[n*n];
		int[] triangles	= new int[(n-1)*(n-1)*6];
		for(int j=0; j<n; j++)
		for(int i=0; i<n; i++)
		{
			X[j*n+i] =new Vector3(5-10.0f*i/(n-1), 0, 5-10.0f*j/(n-1));
			UV[j*n+i]=new Vector3(i/(n-1.0f), j/(n-1.0f));
		}
        int t=0;
		for(int j=0; j<n-1; j++)
		for(int i=0; i<n-1; i++)	
		{
			triangles[t*6+0]=j*n+i;
			triangles[t*6+1]=j*n+i+1;
			triangles[t*6+2]=(j+1)*n+i+1;
			triangles[t*6+3]=j*n+i;
			triangles[t*6+4]=(j+1)*n+i+1;
			triangles[t*6+5]=(j+1)*n+i;
			t++;
		}
		mesh.vertices=X;
		mesh.triangles=triangles;
		mesh.uv = UV;
		mesh.RecalculateNormals ();


		//Construct the original E
		int[] _E = new int[triangles.Length*2];
		for (int i=0; i<triangles.Length; i+=3) 
		{
			_E[i*2+0]=triangles[i+0];
			_E[i*2+1]=triangles[i+1];
			_E[i*2+2]=triangles[i+1];
			_E[i*2+3]=triangles[i+2];
			_E[i*2+4]=triangles[i+2];
			_E[i*2+5]=triangles[i+0];
		}
		//Reorder the original edge list
		for (int i=0; i<_E.Length; i+=2)
			if(_E[i] > _E[i + 1]) 
				Swap(ref _E[i], ref _E[i+1]);
		//Sort the original edge list using quicksort
		Quick_Sort (ref _E, 0, _E.Length/2-1);

		int e_number = 0;
		for (int i=0; i<_E.Length; i+=2)
			if (i == 0 || _E [i + 0] != _E [i - 2] || _E [i + 1] != _E [i - 1]) 
					e_number++;

		E = new int[e_number * 2];
		for (int i=0, e=0; i<_E.Length; i+=2)
			if (i == 0 || _E [i + 0] != _E [i - 2] || _E [i + 1] != _E [i - 1]) 
			{
				E[e*2+0]=_E [i + 0];
				E[e*2+1]=_E [i + 1];
				e++;
            }

		L = new float[E.Length/2];
		for (int e=0; e<E.Length/2; e++) 
		{
			int v0 = E[e*2+0];
			int v1 = E[e*2+1];
			L[e]=(X[v0]-X[v1]).magnitude;
		}

		V = new Vector3[X.Length];
		for (int i=0; i<V.Length; i++)
			V[i] = new Vector3 (0, 0, 0);
    }

    void Quick_Sort(ref int[] a, int l, int r)
	{
		int j;
		if(l<r)
		{
			j=Quick_Sort_Partition(ref a, l, r);
			Quick_Sort (ref a, l, j-1);
			Quick_Sort (ref a, j+1, r);
		}
	}

	int  Quick_Sort_Partition(ref int[] a, int l, int r)
	{
		int pivot_0, pivot_1, i, j;
		pivot_0 = a [l * 2 + 0];
		pivot_1 = a [l * 2 + 1];
		i = l;
		j = r + 1;
		while (true) 
		{
			do ++i; while( i<=r && (a[i*2]<pivot_0 || a[i*2]==pivot_0 && a[i*2+1]<=pivot_1));
			do --j; while(  a[j*2]>pivot_0 || a[j*2]==pivot_0 && a[j*2+1]> pivot_1);
			if(i>=j)	break;
			Swap(ref a[i*2], ref a[j*2]);
			Swap(ref a[i*2+1], ref a[j*2+1]);
		}
		Swap (ref a [l * 2 + 0], ref a [j * 2 + 0]);
		Swap (ref a [l * 2 + 1], ref a [j * 2 + 1]);
		return j;
	}

	void Swap(ref int a, ref int b)
	{
		int temp = a;
		a = b;
		b = temp;
	}

	void Collision_Handling()
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] X = mesh.vertices;

        float radius = 2.7f;
        float tinv = 1.0f / t;
        Vector3 c = sphere.transform.position;
        //Handle colllision.
        for (int i = 0; i < X.Length; ++i){
            float dist = (X[i] - c).magnitude;
			if(dist < radius){
                Vector3 factor = radius * Vector3.Normalize(X[i] - c);
                V[i] += tinv * (c + factor - X[i]);
				X[i] = c + factor;
            }
        }

        mesh.vertices = X;
    }


    
	void Get_Gradient(Vector3[] X, Vector3[] X_hat, float t, Vector3[] G)
	{


        float t2inv = 1.0f / (t * t);
        Vector3 gravity = new Vector3(0.0f, 9.8f, 0.0f) * mass;
        //Momentum and Gravity.
        for (int i = 0; i < X.Length; ++i){
            G[i] = t2inv * mass * (X[i] - X_hat[i]);
            G[i] += gravity;
        }

        //Spring Force.
        for (int i = 0; i < E.Length; i += 2)
        {
            int j = i + 1;
            Vector3 deltaX = X[E[i]] - X[E[j]];
            Vector3 factor = spring_k * (1 - Mathf.Sqrt(L[i / 2] / deltaX.magnitude)) * deltaX;
            G[E[i]] += factor;
            G[E[j]] -= factor;
        }
        G[0] = Vector3.zero;
        G[n - 1] = Vector3.zero;

    }
    bool init = false;
    Vector3[] X_hat, G;
    // Update is called once per frame
    void Update () 
	{
		
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] X 		= mesh.vertices;
        Vector3[] last_X 	= new Vector3[X.Length];
        if (!init)
        {
            X_hat = new Vector3[X.Length];
            G = new Vector3[X.Length];
            init = true;
        }

        //Initial Setup.
        for (int i = 0; i < V.Length; ++i)
        {
            V[i] *= damping;
			// x[0] + t * v[0]
            X_hat[i] = X[i] + t * V[i];
            X[i] = X_hat[i];
        }
        float t2inv = 1.0f / (t * t);
        float omega = 0;
		
        // bool end = false;
        for(int k=0; k<32; k++)
		{
			Get_Gradient(X, X_hat, t, G);
            float ggradinv = 1.0f / (t2inv * mass + 4 * spring_k);

            float ggrad = t2inv * mass + 4 * spring_k;

            switch(k){
				case 0:omega = 1.0f;break;
				case 1:omega = 2.0f / (2.0f - rho * rho);break;
				default:omega = 4.0f / (4.0f - rho * rho * omega);break;
            }
            //Update X by gradient.
            for (int i = 0; i < X.Length; ++i){
				if(i == 0 || i == n - 1) continue;
                Vector3 deltaX = -ggradinv * G[i];

                // Chebyshev Acceleration
				// if rho is closer 1, the algorithm will lose stability.
                // Vector3 r = -G[i] - ggrad * deltaX;
				// if(r.magnitude < 0.1f){
                //     end = true;
                // }
                // Vector3 old_deltaX = deltaX;
                // deltaX = deltaX + (ggradinv * r);
                // deltaX = omega * deltaX + (1 - omega) * last_X[i];
                // last_X[i] = old_deltaX;

                X[i] = X[i] + deltaX;
            }

			// if(end){
            //     break;
            // }
        }

		//Finishing.
		
		mesh.vertices = X;
        float tinv = 1.0f / t;
        for (int i = 0; i < V.Length; ++i){
            V[i] += (X[i] - X_hat[i]) * tinv;
        }
        V[0] = V[n - 1] = Vector3.zero;

        Collision_Handling();
        mesh.RecalculateNormals();
	}
}
