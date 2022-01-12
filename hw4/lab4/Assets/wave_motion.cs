using UnityEngine;
using System.Collections;

public class wave_motion : MonoBehaviour 
{
	int size 		= 100;
	float rate 		= 0.005f;
	float gamma		= 0.001f;
	float damping 	= 0.996f;

	float[,] 	old_h;
	float[,]	low_h;
	float[,]	vh;
	float[,]	b;

	bool [,]	cg_mask;
	float[,]	cg_p;
	float[,]	cg_r;
	float[,]	cg_Ap;
	bool 	tag=true;

	float dt = 0.015f;

	Vector3 	cube_v = Vector3.zero;
	Vector3 	cube_w = Vector3.zero;

    float[,] water_force;
    float dA = 0.01f;
    float gravity = 9.8f;

    // ---------------------------------- block fields ------------------------------------
    float v_damping = 0.999f;
    Transform[] block;
    Vector3[] v = new Vector3[] { new Vector3(0, 0, 0), new Vector3(0, 0, 0) };
    Vector3[] w = new Vector3[] { new Vector3(0, 0, 0), new Vector3(0, 0, 0) };

	Vector3[] force = new Vector3[] { new Vector3(0, 0, 0), new Vector3(0, 0, 0) };

    Matrix4x4[] I_ref = new Matrix4x4[] { Matrix4x4.zero, Matrix4x4.zero };							// reference inertia

    [Range(1.0f, 1000.0f)]
    public float block_mass = 500.0f;
// -----------------------------------------------------------------------------------

    // Use this for initialization
    void Start () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		mesh.Clear ();

		Vector3[] X=new Vector3[size*size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			X[i*size+j].x=i*0.1f-size*0.05f;
			X[i*size+j].y=0;
			X[i*size+j].z=j*0.1f-size*0.05f;
		}

		int[] T = new int[(size - 1) * (size - 1) * 6];
		int index = 0;
		for (int i=0; i<size-1; i++)
		for (int j=0; j<size-1; j++)
		{
			T[index*6+0]=(i+0)*size+(j+0);
			T[index*6+1]=(i+0)*size+(j+1);
			T[index*6+2]=(i+1)*size+(j+1);
			T[index*6+3]=(i+0)*size+(j+0);
			T[index*6+4]=(i+1)*size+(j+1);
			T[index*6+5]=(i+1)*size+(j+0);
			index++;
		}
		mesh.vertices  = X;
		mesh.triangles = T;
		mesh.RecalculateNormals ();

		low_h 	= new float[size,size];
		old_h 	= new float[size,size];
		vh 	  	= new float[size,size];
		b 	  	= new float[size,size];

		cg_mask	= new bool [size,size];
		cg_p 	= new float[size,size];
		cg_r 	= new float[size,size];
		cg_Ap 	= new float[size,size];

        water_force = new float[size, size];

		for (int i=0; i<size; i++)
		for (int j=0; j<size; j++) 
		{
			low_h[i,j]=99999;
			old_h[i,j]=0;
			vh[i,j]=0;
		}

        // if(GameObject.Find("Block") == null || GameObject.Find("Cube") == null){
        //     Debug.LogError("Block Missing");
        // }else{
        //     block1 = GameObject.Find("Block").transform;
        //     block2 = GameObject.Find("Cube").transform;
        // }

		// initialize rigid

        block = new Transform[2];
        block[0] = GameObject.Find("Block").transform;
        block[1] = GameObject.Find("Cube").transform;
    }

	void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
	{
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			Ax[i,j]=0;
			if(i!=0)		Ax[i,j]-=x[i-1,j]-x[i,j];
			if(i!=size-1)	Ax[i,j]-=x[i+1,j]-x[i,j];
			if(j!=0)		Ax[i,j]-=x[i,j-1]-x[i,j];
			if(j!=size-1)	Ax[i,j]-=x[i,j+1]-x[i,j];
		}
	}

	float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
	{
		float ret=0;
		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			ret+=x[i,j]*y[i,j];
		}
		return ret;
	}

	void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
	{
		//Solve the Laplacian problem by CG.
		A_Times(mask, x, cg_r, li, ui, lj, uj);

		for(int i=li; i<=ui; i++)
		for(int j=lj; j<=uj; j++)
		if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
		{
			cg_p[i,j]=cg_r[i,j]=b[i,j]-cg_r[i,j];
		}

		float rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);

		for(int k=0; k<128; k++)
		{
			if(rk_norm<1e-10f)	break;
			A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
			float alpha=rk_norm/Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				x[i,j]   +=alpha*cg_p[i,j];
				cg_r[i,j]-=alpha*cg_Ap[i,j];
			}

			float _rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);
			float beta=_rk_norm/rk_norm;
			rk_norm=_rk_norm;

			for(int i=li; i<=ui; i++)
			for(int j=lj; j<=uj; j++)
			if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
			{
				cg_p[i,j]=cg_r[i,j]+beta*cg_p[i,j];
			}
		}

	}

    void processBlockToWater(int ind, float[,] h, float[,] new_h){
        Mesh blockMesh = block[ind].GetComponent<MeshFilter>().mesh;
        Vector3[] planeX = GetComponent<MeshFilter>().mesh.vertices;
        Bounds aabb = new Bounds(block[ind].position, blockMesh.bounds.size);
        int minX = Mathf.Max((int)(aabb.min.x * 10.0f + size * 0.5f), 0);
        int maxX = Mathf.Min((int)(aabb.max.x * 10.0f + size * 0.5f), size - 1);
        int minY = Mathf.Max((int)(aabb.min.z * 10.0f + size * 0.5f), 0);
        int maxY = Mathf.Min((int)(aabb.max.z * 10.0f + size * 0.5f), size - 1);
        // Debug.Log("minX: " + minX + ", maxX: " + maxX);
        // Debug.Log("minY: " + minY + ", maxY: " + maxY);
        float e = Mathf.Abs(block[ind].position.y - 0.5f);
        for (int i = minX; i <= maxX;++i)
        {
            for (int j = minY; j <= maxY;++j) {
                Vector3 origin = planeX[i * size + j];
                origin.y = transform.position.y;
                Vector3 dir = block[ind].position - origin;
                if (!Physics.Raycast(origin, dir, 1000.0f)) {
                    low_h[i, j] = h[i, j] - e;
                   	b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
                    cg_mask[i, j] = true;
                } else {
                    cg_mask[i, j] = false;
                }
            }
        }

        Conjugate_Gradient(cg_mask, b, vh, minX, maxX, minY, maxY);

    }

    void processWaterToBlock(int ind, float[,] h, float[,] new_h){
        Mesh blockMesh = block[ind].GetComponent<MeshFilter>().mesh;
        Bounds aabb = new Bounds(block[ind].position, blockMesh.bounds.size);
        int minX = Mathf.Max((int)(aabb.min.x * 10.0f + size * 0.5f), 0);
        int maxX = Mathf.Min((int)(aabb.max.x * 10.0f + size * 0.5f), size - 1);
        int minY = Mathf.Max((int)(aabb.min.z * 10.0f + size * 0.5f), 0);
        int maxY = Mathf.Min((int)(aabb.max.z * 10.0f + size * 0.5f), size - 1);
        float e = Mathf.Abs(block[ind].position.y - 0.5f);
		for (int i = minX; i <= maxX;++i)
        {
            for (int j = minY; j <= maxY;++j) {
                if(cg_mask[i,j]){
                    water_force[i, j] = 1000 * gravity * dA * e;
                }else{
                    water_force[i, j] = 0;
                }
            }
        }
	}

	void Shallow_Wave(float[,] old_h, float[,] h, float [,] new_h)
	{
        //Step 1:
        //TODO: Compute new_h based on the shallow wave model.
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                new_h[i, j] = h[i, j] + (h[i, j] - old_h[i, j]) * damping
                    + (
						(i == 0?h[i, j] :h[i - 1, j]) 
						+ (i == size - 1 ? h[i, j] : h[i + 1, j])
                        + (j == 0 ? h[i, j] : h[i, j - 1])
                        + (j == size - 1 ? h[i, j] : h[i, j + 1])
						- 4 * h[i, j]
						) * rate;
            }
        }

        //Step 2: Block->Water coupling
        //TODO: for block 1, calculate low_h.
        //TODO: then set up b and cg_mask for conjugate gradient.
        //TODO: Solve the Poisson equation to obtain vh (virtual height).
        processBlockToWater(0, h, new_h);

        //TODO: for block 2, calculate low_h.
        //TODO: then set up b and cg_mask for conjugate gradient.
        //TODO: Solve the Poisson equation to obtain vh (virtual height).
        processBlockToWater(1, h, new_h);

        //TODO: Diminish vh.
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                vh[i, j] = gamma * vh[i, j];
            }
        }
        //TODO: Update new_h by vh.
        for (int i = 0; i < size; ++i)
        {
            for (int j = 0; j < size; ++j)
            {
                new_h[i, j] += (
                    (i == 0 ? vh[i, j] : vh[i - 1, j])
                    + (i == size - 1 ? vh[i, j] : vh[i + 1, j])
                    + (j == 0 ? vh[i, j] : vh[i, j - 1])
                    + (j == size - 1 ? vh[i, j] : vh[i, j + 1])
                    - 4 * vh[i, j]
					) * rate;
            }
        }

        // for (int i = 0; i < size; ++i)
        // {
        //     for (int j = 0; j < size; ++j)
        //     {
        //         if(cg_mask[i,j]){
        //             // Debug.Log("vh: " + vh[i, j]);
        //             water_force[i, j] = 1000 * gravity * dA * (vh[i,j] * gamma);
        //             // Debug.Log("water_force:" + water_force[i, j]);
        //         }else{
        //             water_force[i, j] = 0;
        //         }
        //     }
        // }


        //Step 3
        //TODO: old_h <- h; h <- new_h;
        for (int i = 0; i < size;++i){
            for (int j = 0; j < size; ++j)
            {
                old_h[i, j] = h[i, j];
                h[i, j] = new_h[i, j];
                vh[i, j] = 0;
            }
        }

        //Step 4: Water->Block coupling.
        //More TODO here.
        processWaterToBlock(0, h, new_h);
        processWaterToBlock(1, h, new_h);
    }

    const float EPS = 1e-4f;
    void Update () 
	{
		Mesh mesh = GetComponent<MeshFilter> ().mesh;
		Vector3[] X    = mesh.vertices;
		float[,] new_h = new float[size, size];
		float[,] h     = new float[size, size];

        //TODO: Load X.y into h.
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                h[i, j] = X[i * size + j].y;
            }
        }

        if (Input.GetKeyDown ("r")) 
		{
            //TODO: Add random water.
            float height = Random.Range(0.01f, 0.5f);
            int x = Random.Range(0, size - 1);
            int y = Random.Range(0, size - 1);
            h[x, y] += height;
			if(x > 0 && x < size - 1 && y > 0 && y < size - 1){
                float r = height * 0.25f;
                h[x - 1, y] -= r;
                h[x + 1, y] -= r;
                h[x, y - 1] -= r;
                h[x, y + 1] -= r;
			}else if(x == 0){
                float r = height * 0.5f;
				if(y == 0){
                    h[x + 1, y] -= r;
                    h[x, y + 1] -= r;
				}else{
                    h[x + 1, y] -= r;
                    h[x, y - 1] -= r;
				}
			}else{
                float r = height * 0.5f;
                if (y == 0){
                    h[x - 1, y] -= r;
                    h[x, y + 1] -= r;
                }
                else{
                    h[x - 1, y] -= r;
                    h[x, y - 1] -= r;
                }
			}
        }

        // update rigid
        for (int i = 0; i < 2; ++i)
        {
            Debug.Log(force[i]);
            Vector3 v1;
            v1 = v[i] + dt * force[i] / block_mass;
            block[i].position = block[i].position + dt * v1;

            v[i] = v1 * damping;
            Debug.Log("v: " + v[i]);
			if(v[i].magnitude < EPS){
                v[i] = Vector3.zero;
            }
            force[i] = Vector3.zero;
        }

		for(int l=0; l<8; l++)
		{
			Shallow_Wave(old_h, h, new_h);
        }

        for (int ind = 0; ind < 2; ++ind)
        {
            Mesh blockMesh = block[ind].GetComponent<MeshFilter>().mesh;
            Vector3[] planeX = GetComponent<MeshFilter>().mesh.vertices;
            Bounds aabb = new Bounds(block[ind].position, blockMesh.bounds.size);
            int minX = Mathf.Max((int)(aabb.min.x * 10.0f + size * 0.5f), 0);
            int maxX = Mathf.Min((int)(aabb.max.x * 10.0f + size * 0.5f), size - 1);
            int minY = Mathf.Max((int)(aabb.min.z * 10.0f + size * 0.5f), 0);
            int maxY = Mathf.Min((int)(aabb.max.z * 10.0f + size * 0.5f), size - 1);
            for (int i = minX; i <= maxX; ++i)
            {
                for (int j = minY; j <= maxY; ++j)
                {
                    if (cg_mask[i, j])
                    {
                        // Debug.Log("water_force:" + water_force[i, j]);
                        force[ind].y += water_force[i, j];
                    }
                    else
                    {
                        water_force[i, j] = 0;
                    }
                }
            }
            force[ind].y -= gravity * block_mass;
        }
        cg_mask = new bool[size, size];

        //TODO: Store h back into X.y and recalculate normal.

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < size; j++)
            {
                X[i * size + j].y = h[i, j];
            }
        }
        mesh.vertices = X;
        mesh.RecalculateNormals();
		

    }
}
