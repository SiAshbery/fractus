using System.Collections;
using System.Collections.Generic;
using UnityEngine;

[RequireComponent(typeof(Camera))]
[ExecuteInEditMode]
public class RaymarchCamera : SceneViewFilter
{
    // Create raymarch screenspace shader.
    [SerializeField]
    private Shader _shader;

    public Material _raymarchMaterial
    {
        get
        {
            if (!_raymarchMat && _shader)
            {
                _raymarchMat = new Material(_shader);
                _raymarchMat.hideFlags = HideFlags.HideAndDontSave;
            }
            return _raymarchMat;
        }
    }

    private Material _raymarchMat;

    // Get raymarch camera.
    public Camera _camera
    {
        get
        {
            if (!_cam)
            {
                _cam = GetComponent<Camera>();
            }
            return _cam;
        }
    }

    private Camera _cam;

    // Use editor directional light to calculate shader lighting
    [Header("Lighting")]
    public Light _directionLight;
    public Color _mainCol;
    public float _lightIntensity, _shadowIntensity, _shadowPenumbra;
    public Vector2 _shadowDistance;

    [Header("Ambient Occlusion")]
    public float _aoStepSize;
    public float _aoIntensity;
    public int _aoIterations;

    [Header("Fog")]
    public Color _fogColor;
    [Range(0.0f, 1.0f)]
    public float _fogDensity;

    [Header("Raymarching Controls")]
    // Set the max distance of ray travel in editor
    public float _maxDistance;
    // Set the max iterations of ray travel in editor
    public int _maxIterations;
    // set the acurracy of ray marching hits
    [Range(0.1f, 0.001f)]
    public float _accuracy;

    // Mod repetition interval
    public Vector3 _modInterval;

    [Header("Sphere")]
    public Vector4 _sphere1;

    [Header("Box")]
    public Vector4 _box1;

    [Header("Sierpinski Triangle")]
    public Vector3 _recursiveTet1;
    public float _recursiveTet1Offset;
    public int _recursiveTet1Iterations;

    [Header("Mandel Bulb")]
    public Vector3 _mandelBulb1;
    public float _mandelBulb1Power, _mandelBulb1Bailout;
    public int _mandelBulb1Iterations;

    [Header("Mandel Box")]
    public Vector3 _mandelBox1;
    public Vector3 _mandelBox1FoldLimit;
    public float _mandelBox1Scale, _mandelBox1SphereRadius;
    public int _mandelBox1Iterations;

    [Header("Apollonian")]
    public Vector3 _apollonian1;
    public float _apollonian1Scale;
    public int _apollonian1Iterations;
    public Vector3 _apollonian1Size;

    // On render image to communicate with shader
    private void OnRenderImage(RenderTexture source, RenderTexture destination)
    {
        // check if raymarch metarial has been set.
        // if it hasn't just return the image
        if (!_raymarchMaterial)
        {
            // Copies source texture into destination render texture with a shader.
            //This is mostly used for implementing post-processing effects.
            Graphics.Blit(source, destination);
            return;
        }

        //Set variables to shader
        _raymarchMaterial.SetMatrix("_CamFrustum", CamFrustum(_camera));
        _raymarchMaterial.SetMatrix("_CamToWorld", _camera.cameraToWorldMatrix);
        _raymarchMaterial.SetFloat("_maxDistance", _maxDistance);
        _raymarchMaterial.SetInt("_maxIterations", _maxIterations);
        _raymarchMaterial.SetFloat("_accuracy", _accuracy);
        // Lighting
        _raymarchMaterial.SetVector("_lightDir", _directionLight ? _directionLight.transform.forward : Vector3.down);
        _raymarchMaterial.SetColor("_lightCol", _directionLight ?_directionLight.color : new Color(1,1,1));
        _raymarchMaterial.SetColor("_mainCol", _mainCol);
        _raymarchMaterial.SetColor("_fogColor", _fogColor);
        _raymarchMaterial.SetFloat("_fogDensity", _fogDensity);
        _raymarchMaterial.SetFloat("_lightIntensity", _lightIntensity);
        _raymarchMaterial.SetFloat("_shadowIntensity", _shadowIntensity);
        _raymarchMaterial.SetVector("_shadowDistance", _shadowDistance);
        _raymarchMaterial.SetFloat("_shadowPenumbra", _shadowPenumbra);
        // AO
        _raymarchMaterial.SetFloat("_aoStepSize", _aoStepSize);
        _raymarchMaterial.SetFloat("_aoIntensity", _aoIntensity);
        _raymarchMaterial.SetInt("_aoIterations", _aoIterations);
        // Shapes
        _raymarchMaterial.SetVector("_sphere1", _sphere1);
        _raymarchMaterial.SetVector("_box1", _box1);
        // Recursive tetrahedron
        _raymarchMaterial.SetVector("_recursiveTet1", _recursiveTet1);
        _raymarchMaterial.SetFloat("_recursiveTet1Offset", _recursiveTet1Offset);
        _raymarchMaterial.SetInt("_recursiveTet1Iterations", _recursiveTet1Iterations);
        // Mandel Bulb
        _raymarchMaterial.SetVector("_mandelBulb1", _mandelBulb1);
        _raymarchMaterial.SetInt("_mandelBulb1Iterations", _mandelBulb1Iterations);
        _raymarchMaterial.SetFloat("_mandelBulb1Power", _mandelBulb1Power);
        _raymarchMaterial.SetFloat("_mandelBulb1Bailout", _mandelBulb1Bailout);
        // Mandel Box
        _raymarchMaterial.SetVector("_mandelBox1", _mandelBox1);
        _raymarchMaterial.SetInt("_mandelBox1Iterations", _mandelBox1Iterations);
        _raymarchMaterial.SetFloat("_mandelBox1Scale", _mandelBox1Scale);
        _raymarchMaterial.SetVector("_mandelBox1FoldLimit", _mandelBox1FoldLimit);
        _raymarchMaterial.SetFloat("_mandelBox1SphereRadius", _mandelBox1SphereRadius);
        // Apollonian
        _raymarchMaterial.SetVector("_apollonian1", _apollonian1);
        _raymarchMaterial.SetFloat("_apollonian1Scale", _apollonian1Scale);
        _raymarchMaterial.SetInt("_apollonian1Iterations", _apollonian1Iterations);
        _raymarchMaterial.SetVector("_apollonian1Size", _apollonian1Size);
    // Repetition
    _raymarchMaterial.SetVector("_modInterval", _modInterval);

        // Here we create the quad which represents our screen.
        // It is important that it's corner corresponds correctly to the corners of the frustum.
        RenderTexture.active = destination;

        // Show skybox and standard unity scene meshes with raymarched shapes.
        // To do this we assign the main tex of the shader to the current source
        // Then we can apply the ray marching shader and return a combination of both.
        _raymarchMaterial.SetTexture("_MainTex", source);

        // Push as quad into the text coord 0 of the shader.
        GL.PushMatrix();
        // Load the orthographic view
        GL.LoadOrtho();
        // Set pass of the raymarch material to 0
        _raymarchMaterial.SetPass(0);
        // Start drawing the quad
        GL.Begin(GL.QUADS);

        // We will start at the bottom left and move to each corner in sequence
        // BL
        GL.MultiTexCoord2(0, 0.0f, 0.0f);
        // draw a vertex.
        // NOTE the third argument is not needed as we are in orthographic view,
        // but we can use it to align the vertex with the corresponding corner of the frustum
        GL.Vertex3(0.0f, 0.0f, 3.0f);
        // BR
        GL.MultiTexCoord2(0, 1.0f, 0.0f);
        GL.Vertex3(1.0f, 0.0f, 2.0f);
        // TR
        GL.MultiTexCoord2(0, 1.0f, 1.0f);
        GL.Vertex3(1.0f, 1.0f, 1.0f);
        // TL
        GL.MultiTexCoord2(0, 0.0f, 1.0f);
        GL.Vertex3(0.0f, 1.0f, 0.0f);

        GL.End();
        GL.PopMatrix();
    }

    //Camera Frustum
    // There are two things we need toset:
    // - Camera frustum which consists of 4 directions that correspond to the four corners of the screen.
    // - A render texture quad on which we will render the output of the rays as a render texture.
    //We need to make sure that all corners of the quad correspond to the four directions of the frustum.

    // to create the frustum matrix we need to get camer's field-of-view and it's aspect ratio
    // With these values we can calculate each corner

    // The first thing we need to do is get the tangent of the field of view in radians divided by 2
    // We will then create a vector to go up (goUp) and one to go right (goRight).
    // We then calculate the four courners of the frustum starting with the top left corner
    // First we move backwards which is equal to -Vector3.forward and then we go left which is equal to -goRight and then we add goUp
    // We complete the other corners by following the same steps adjust each movement accordingly and we store themin a matrix.

    private Matrix4x4 CamFrustum(Camera cam)
    {
        Matrix4x4 frustum = Matrix4x4.identity;
        // get tangent of FOV in radians
        float fov = Mathf.Tan((cam.fieldOfView * 0.5f) * Mathf.Deg2Rad);

        // To get all corners of the frustum we create two vectors
        // One to go up and one to go right
        Vector3 goUp = Vector3.up * fov;
        Vector3 goRight = Vector3.right * fov * cam.aspect;

        //TOP LEFT
        Vector3 TL = (-Vector3.forward - goRight + goUp);
        //TOP RIGHT
        Vector3 TR = (-Vector3.forward + goRight + goUp);
        //BOTTOM RIGHT
        Vector3 BR = (-Vector3.forward + goRight - goUp);
        //BOTTOM LEFT
        Vector3 BL = (-Vector3.forward - goRight - goUp);

        // Apply vectors to frustum matrix
        frustum.SetRow(0, TL);
        frustum.SetRow(1, TR);
        frustum.SetRow(2, BR);
        frustum.SetRow(3, BL);

        return frustum;

        // Now in the on render image we can set the matrix of the shader to the frustum matrix
    }

}
