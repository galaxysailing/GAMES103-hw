using UnityEngine;
using System.Collections;

public class Rigid_Bunny : MonoBehaviour 
{
	bool launched 		= false;
	float dt 			= 0.015f;
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity
	
	float mass;									// mass
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.999f;				// for velocity decay
	float angular_decay	= 0.98f;				
	float restitution 	= 0.5f;                 // for collision

	// gravity force
    Vector3 g_force = new Vector3(0.0f, 0.0f, 0.0f);

    Vector3[] vertex_force;


    ///////////////////// Math helper //////////////////////////

    Quaternion add(Quaternion a, Quaternion b)
    {
        return Quaternion.Normalize(new Quaternion(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w));
    }
	
	Matrix4x4 subtract(Matrix4x4 a, Matrix4x4 b){
        Matrix4x4 res = Matrix4x4.zero;
        res[0, 0] = a[0, 0] - b[0, 0];
        res[0, 1] = a[0, 1] - b[0, 1];
        res[0, 2] = a[0, 2] - b[0, 2];
        res[0, 3] = a[0, 3] - b[0, 3];

        res[1, 0] = a[1, 0] - b[1, 0];
        res[1, 1] = a[1, 1] - b[1, 1];
        res[1, 2] = a[1, 2] - b[1, 2];
        res[1, 3] = a[1, 3] - b[1, 3];

        res[2, 0] = a[2, 0] - b[2, 0];
        res[2, 1] = a[2, 1] - b[2, 1];
        res[2, 2] = a[2, 2] - b[2, 2];
        res[2, 3] = a[2, 3] - b[2, 3];

        res[3, 0] = a[3, 0] - b[3, 0];
        res[3, 1] = a[3, 1] - b[3, 1];
        res[3, 2] = a[3, 2] - b[3, 2];
        res[3, 3] = a[3, 3] - b[3, 3];

        return res;
    }

    //////////////////////////////////////////////////////////////////
    

    // Use this for initialization
    void Start () 
	{		
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;
        vertex_force = new Vector3[vertices.Length];
		
        // gravity accelerated speed
        Vector3 g = new Vector3(0.0f, -9.8f, 0.0f);

        float m=1;
		mass=0;
		for (int i=0; i<vertices.Length; i++) 
		{
			mass += m;
			float diag=m*vertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
			I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
			I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
			I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
			I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
			I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
			I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
			I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
			I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];

            vertex_force[i] = Vector3.zero;
        }
        I_ref[3, 3] = 1;
        g_force = mass * g;
    }
	
	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0; 
		A [0, 1] = -a [2]; 
		A [0, 2] = a [1]; 
		A [1, 0] = a [2]; 
		A [1, 1] = 0; 
		A [1, 2] = -a [0]; 
		A [2, 0] = -a [1]; 
		A [2, 1] = a [0]; 
		A [2, 2] = 0; 
		A [3, 3] = 1;
		return A;
	}

    //////////////////// Collision Helper ////////////////////

    const float EPS = 0.05f;
    bool Detect_Collision(Vector3 P, Vector3 N, out Vector3 vi, out Vector3 ri){
        Matrix4x4 R = Matrix4x4.Rotate(transform.rotation);
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] vertices = mesh.vertices;
        ri = Vector3.zero;
        vi = Vector3.zero;
        int cnt = 0;
        for (int i = 0; i < vertices.Length; ++i)
        {
            Vector3 pos = new Vector4(transform.position.x, transform.position.y, transform.position.z, 0.0f) + (R * vertices[i]);
            float phi = Vector3.Dot(pos - P, N);
            if (phi <= 0)
            {
                Vector3 tmp = v + Vector3.Cross(w, R * vertices[i]);
                if(Vector3.Dot(tmp, N) < 0){
                    vi += tmp - v;
                    ri += vertices[i];
                    ++cnt;
                }
            }
        }
		if(cnt == 0){
            return false;
        }
        vi /= cnt;
        vi += v;
        ri /= cnt;
        return true;
    }

    //////////////////////////////////////////////////////////
    // In this function, update v and w by the impulse due to the collision with
    //a plane <P, N>
    void Collision_Impulse(Vector3 P, Vector3 N)
	{
        Vector3 vi, ri;
        if(!Detect_Collision(P, N, out vi, out ri)){
            return;
        }

        // compute
        Matrix4x4 R = Matrix4x4.Rotate(transform.rotation);
        Matrix4x4 invI = (R * I_ref * R.transpose).inverse;
        Vector3 vn = Vector3.Dot(vi, N) * N;
        Vector3 vt = vi - vn;
        float a = Mathf.Max(0, 1 - restitution * (1 + restitution) * vn.sqrMagnitude / vt.sqrMagnitude);
        vn = -restitution * vn;
        vt = a * vt;
        Vector3 vi_new = vn + vt;
        Matrix4x4 Rri = Get_Cross_Matrix(R * ri);
        Matrix4x4 Kprev = Matrix4x4.identity;

        float inv_mass = 1.0f / mass;
        Kprev[0, 0] *= inv_mass;
        Kprev[1, 1] *= inv_mass;
        Kprev[2, 2] *= inv_mass;
        Kprev[3, 3] *= inv_mass;
        Matrix4x4 invK = subtract(Kprev, Rri * invI * Rri).inverse;

        Vector3 j = invK * (vi_new - vi);

        v = v + j * inv_mass;
        w = new Vector4(w.x, w.y, w.z, 0.0f) + invI * (Rri * j);
    }

	void Update_Velocities(out Vector3 v_mid){
		
		// update linear velocity
        Vector3 v0 = v;
        if (v0.magnitude < EPS)
        {
            v = Vector3.zero;
            v_mid = Vector3.zero;
        } else {
            v = v0 + dt * (g_force / mass);
            v = linear_decay * v;
            v_mid = (v + v0) * 0.5f;
        }

        // update angular velocity
        Vector3 w0 = w;

        if (w0.magnitude < EPS) {
            w = Vector3.zero;
        } else {
            Matrix4x4 R = Matrix4x4.Rotate(transform.rotation);
            Matrix4x4 I = R * I_ref * R.transpose;
            Vector3 torque = Vector3.zero;
            Mesh mesh = GetComponent<MeshFilter>().mesh;
            Vector3[] vertices = mesh.vertices;
            for (int i = 0; i < vertices.Length; ++i)
            {
                torque += Vector3.Cross(R * vertices[i], new Vector3(0.0f, -9.8f, 0.0f));
            }
            Matrix4x4 inv_I = I.inverse;
            w = new Vector4(w0.x, w0.y, w0.z, 0.0f) + inv_I * torque * dt;
            w = angular_decay * w;
        }

    }

    // Update is called once per frame
    void Update () 
	{
        //Game Control
        if(Input.GetKey("r"))
		{
			transform.position = new Vector3(0, 0.6f, 0);
			restitution = 0.5f;
			launched=false;
		}
		if(Input.GetKey("l"))
		{
			v = new Vector3 (5, 2, 0);
            // w = new Vector3(50, 20, 0);
            launched=true;
		}

		if(!launched){
            return;
        }

        Vector3 v_mid = v;
        // Part I: Update velocities
        Update_Velocities(out v_mid);

        // Part II: Collision Impulse
        Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
        Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));

        // Part III: Update position & orientation
        //Update linear status
        Vector3 x    = transform.position;
		//Update angular status
		Quaternion q = transform.rotation;

        //leapfrog
        x = x + dt * v_mid;

        q = add(q, new Quaternion(dt * 0.5f * w.x, dt * 0.5f * w.y, dt * 0.5f * w.z, 0.0f) * q);
        // Part IV: Assign to the object
        transform.position = x;
		transform.rotation = q;
	}
}

