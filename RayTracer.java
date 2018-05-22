//package RayTracer;
package RayTracer;

import java.awt.Transparency;
import java.awt.color.*;
import java.awt.image.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Random;
import javax.imageio.ImageIO;


/**
 *  Main class for ray tracing exercise.
 */


public class RayTracer {
	static double epsilon = 0.001;
	static double trigEpsilon = 0.0000001;
	double recMax;// this is now a class variable
	public int imageWidth;
	public int imageHeight;
	public Scene scene;
	
	public enum Shape { SPHERE, PLANE, TRIANGLE , BACKGROUND };
	/**
	 * Runs the ray tracer. Takes scene file, output image file and image size as input.
	 */
	
	
	public static void main(String[] args) {
		
		try {

			RayTracer tracer = new RayTracer();

                        // Default values:
			tracer.imageWidth = 500;
			tracer.imageHeight = 500;

			if (args.length < 2)
				throw new RayTracerException("Not enough arguments provided. Please specify an input scene file and an output image file for rendering.");

			String sceneFileName = args[0];
			String outputFileName = args[1];

			if (args.length > 3)
			{
				tracer.imageWidth = Integer.parseInt(args[2]);
				tracer.imageHeight = Integer.parseInt(args[3]);
			}


			// Parse scene file:
			tracer.scene = tracer.parseScene(sceneFileName);
			if (tracer.scene == null){
				System.out.println("Scene was not properly defined. Program terminates");
			}

			// Render scene:
			tracer.renderScene(outputFileName,tracer);

//		} catch (IOException e) {
//			System.out.println(e.getMessage());
		} catch (RayTracerException e) {
			System.out.println(e.getMessage());
		} catch (Exception e) {
			System.out.println(e.getMessage());
		}


	}

	/**
	 * Parses the scene file and creates the scene. Change this function so it generates the required objects.
	 */
   public Scene parseScene(String sceneFileName) throws IOException, RayTracerException
	{
		FileReader fr = new FileReader(sceneFileName);

		BufferedReader r = new BufferedReader(fr);
		String line = null;
		int lineNum = 0;
		System.out.println("Started parsing scene file " + sceneFileName);
		
		Scene scene= new Scene();
		
		while ((line = r.readLine()) != null)
		{
			line = line.trim();
			++lineNum;

			if (line.isEmpty() || (line.charAt(0) == '#'))
			{  // This line in the scene file is a comment
				continue;
			}
			else
			{
				String code = line.substring(0, 3).toLowerCase();
				// Split according to white space characters:
				String[] params = line.substring(3).trim().toLowerCase().split("\\s+");

				if (code.equals("cam"))
				{
                    // Add code here to parse camera parameters
					int numCamVariables = 11;
					//String[] camString = line.split("\\s+");
					if (params.length != numCamVariables){
						System.out.print("Inadequate number of camera parameters\n");
						return null;
					}
					scene.camera = new Camera(params);
					System.out.println(String.format("Parsed camera parameters (line %d)", lineNum));
				}
				else if (code.equals("set"))
				{
					// Add code here to parse general settings parameters
					int numSetVariables = 5;
					if (params.length < numSetVariables || numSetVariables + 1 < params.length){
						System.out.print("Inadequate number of settings parameters\n");
						return null;
					}
					scene.sceneSettings = new SceneSettings(params);
					recMax = scene.sceneSettings.recMax;
					System.out.println(String.format("Parsed general settings (line %d)", lineNum));
				}
				else if (code.equals("mtl"))
				{
					// Add code here to parse material parameters
					int numMaterialVariables = 11;
					if (params.length != numMaterialVariables){
						System.out.print("Inadequate number of material parameters\n");
						return null;
					}
					// at first material string, add materials array-list.
					if (scene.materials == null)
						scene.materials = new ArrayList<Material>();
					scene.materials.add(new Material(params));
					System.out.println(String.format("Parsed material (line %d)", lineNum));
					
				}
				else if (code.equals("sph"))
				{
					// Add code here to parse sphere parameters
					int numSphereVariables = 5;
					if (params.length != numSphereVariables){
						System.out.print("Inadequate number of sphere parameters\n");
						return null;
					}
					if (scene.spheres == null){
						scene.spheres = new ArrayList<Sphere>();
					}
					scene.spheres.add(new Sphere(params));
					System.out.println(String.format("Parsed sphere (line %d)", lineNum));
				}
				else if (code.equals("pln"))
				{
                    // Add code here to parse plane parameters
					int numPlaneVariables = 5;
					if (params.length != numPlaneVariables){
						System.out.print("Inadequate number of plain parameters\n");
						return null;
					}
					if (scene.planes == null){
						scene.planes = new ArrayList<Plane>();
					}
					scene.planes.add(new Plane(params));
					System.out.println(String.format("Parsed plane (line %d)", lineNum));
				}
				else if (code.equals("trg"))
				{
					// Add code here to parse cylinder parameters
					int numTrigVariables = 10; 
					if (params.length != numTrigVariables){
						System.out.print("Inadequate number of triangle parameters\n");
						return null;
					}
					if (scene.triangles == null){
						scene.triangles = new ArrayList<Triangle>();
					}
					scene.triangles.add(new Triangle(params));
					System.out.println(String.format("Parsed triangle (line %d)", lineNum));
				}
				else if (code.equals("lgt"))
				{
					// Add code here to parse cylinder parameters
					int numLightVariables = 9; 
					if (params.length != numLightVariables){
						System.out.print("Inadequate number of light parameters\n");
						return null;
					}
					if (scene.lights == null){
						scene.lights = new ArrayList<Light>();
					}
					scene.lights.add(new Light(params));
					System.out.println(String.format("Parsed light (line %d)", lineNum));
				}
				else
				{
					System.out.println(String.format("ERROR: Did not recognize object: %s (line %d)", code, lineNum));
				}
			}
		}

        // It is recommended that you check here that the scene is valid,
        // for example camera settings and all necessary materials were defined.
		if (scene.camera == null){
			System.out.println("Invalid scene. Camera parameters were not defined\n");
			return null;
		}
		if (scene.sceneSettings == null){
			System.out.println("Invalid scene. Camera settings were not defined\n");
			return null;
		}
		if (scene.materials == null){
			System.out.println("Invalid scene. No materials were defined\n");
			return null;
		}
		if (scene.lights == null){
			System.out.println("Invalid scene. No lights were defined\n");
			return null;
		}
		if (scene.spheres == null && scene.planes == null && scene.triangles == null){
			System.out.println("Invalid scene. No objectss were defined\n");
			return null;
		}
		System.out.println("Finished parsing scene file " + sceneFileName);
		return scene;
	}

	/**
	 * Renders the loaded scene and saves it to the specified file location.
	 */
	public void renderScene(String outputFileName,RayTracer tracer)
	{
		long startTime = System.currentTimeMillis();
		// Create a byte array to hold the pixel data:
		Random rand = new Random();
		int samplingLevel = tracer.scene.sceneSettings.samplingLevel;
		byte[] rgbData = new byte[this.imageWidth * this.imageHeight * 3];
		for (int i = 0; i < this.imageWidth ; i++){
			for (int j = 0; j < this.imageHeight; j++){
				Vector averagedPixelColor = new Vector();;
				averagedPixelColor.xCor = 0; averagedPixelColor.yCor = 0; averagedPixelColor.zCor = 0;
				for (int Ni = 0; Ni < samplingLevel;Ni++){
					for (int Nj = 0; Nj < samplingLevel;Nj++){
						Ray ray = new Ray(tracer,scene,i,j,Ni,Nj,rand,tracer.imageWidth,tracer.imageHeight);
						IntersectionParameters intersecParam = FindIntersection(ray,tracer.scene,tracer.scene.camera.camPosition,-1);
						Vector PixelColor;
						if (intersecParam.closestShapeType != Shape.BACKGROUND){
							Vector IntesectionPos = getShapePosittion(intersecParam.intersectionDist,ray,tracer.scene.camera.camPosition);
							Vector shapeNormal = getShapeNoraml(tracer.scene,IntesectionPos,intersecParam.closestShapeType,intersecParam.closestShapeIndx);
							Material mat = null;
							PixelColor = getPixelColor(intersecParam.closestShapeType,intersecParam.closestShapeIndx,IntesectionPos,shapeNormal,tracer,rand,mat,ray,0);
						}
						else {
							PixelColor = tracer.scene.sceneSettings.backgroundRGB;
						}
						averagedPixelColor = Vector.addVector(averagedPixelColor, PixelColor);
					}
				}
				averagedPixelColor = Vector.scalarMult(averagedPixelColor, (1/(double)(samplingLevel*samplingLevel)));
				averagedPixelColor = threshhold(averagedPixelColor);
				rgbData[(j * this.imageWidth + i) * 3] = (byte) (averagedPixelColor.xCor*255);
				rgbData[(j * this.imageWidth + i) * 3 + 1] = (byte) (averagedPixelColor.yCor*255);
				rgbData[(j * this.imageWidth + i) * 3 + 2] = (byte) (averagedPixelColor.zCor*255);
			}
		}
		
		long endTime = System.currentTimeMillis();
		Long renderTime = endTime - startTime;

        // The time is measured for your own conveniece, rendering speed will not affect your score
        // unless it is exceptionally slow (more than a couple of minutes)
		System.out.println("Finished rendering scene in " + renderTime.toString() + " milliseconds.");

        // This is already implemented, and should work without adding any code.
		saveImage(this.imageWidth, rgbData, outputFileName);

		System.out.println("Saved file " + outputFileName);

	}


	//////////////////////// FUNCTIONS TO SAVE IMAGES IN PNG FORMAT //////////////////////////////////////////

	/*
	 * Saves RGB data as an image in png format to the specified location.
	 */
	public static void saveImage(int width, byte[] rgbData, String fileName)
	{
		try {

			BufferedImage image = bytes2RGB(width, rgbData);
			ImageIO.write(image, "png", new File(fileName));

		} catch (IOException e) {
			System.out.println("ERROR SAVING FILE: " + e.getMessage());
		}

	}

	/*
	 * Producing a BufferedImage that can be saved as png from a byte array of RGB values.
	 */
	public static BufferedImage bytes2RGB(int width, byte[] buffer) {
	    int height = buffer.length / width / 3;
	    ColorSpace cs = ColorSpace.getInstance(ColorSpace.CS_sRGB);
	    ColorModel cm = new ComponentColorModel(cs, false, false,
	            Transparency.OPAQUE, DataBuffer.TYPE_BYTE);
	    SampleModel sm = cm.createCompatibleSampleModel(width, height);
	    DataBufferByte db = new DataBufferByte(buffer, width * height);
	    WritableRaster raster = Raster.createWritableRaster(sm, db, null);
	    BufferedImage result = new BufferedImage(cm, raster, false, null);

	    return result;
	}
	
	public static IntersectionParameters FindIntersection(Ray ray,Scene scene, Vector raySource,int ShapeIndx){

		
		// create a ray from camera through pixel center
		IntersectionParameters intersecParam = new IntersectionParameters(); 
		intersecParam.intersectionDist = Double.MAX_VALUE;
		// check intersection with all spheres
		if (scene.spheres != null){
			for (int sphereIndx=0; sphereIndx<scene.spheres.size();sphereIndx++){
				double current_t = sphereIntersect(ray,scene.spheres.get(sphereIndx), raySource);
				if (0 <= current_t && current_t < intersecParam.intersectionDist && sphereIndx != ShapeIndx){
					intersecParam.intersectionDist = current_t;
					intersecParam.closestShapeType = Shape.SPHERE;
					intersecParam.closestShapeIndx = sphereIndx;
				}
			}
		}
		// check intersection with all planes
		if (scene.planes != null){
			for (int planeIndex=0; planeIndex<scene.planes.size();planeIndex++){
				double current_t = planeIntersect(ray,scene.planes.get(planeIndex), raySource);
				if (0 <= current_t && current_t < intersecParam.intersectionDist && planeIndex != ShapeIndx){
					intersecParam.intersectionDist = current_t;
					intersecParam.closestShapeType = Shape.PLANE;
					intersecParam.closestShapeIndx = planeIndex;
				}
			}
		}
		// check intersection with all triangles
		if (scene.triangles != null){
			for (int trigIndx=0; trigIndx<scene.triangles.size();trigIndx++){
				double current_t = triangleIntersect(ray,scene.triangles.get(trigIndx),raySource);
				if (Math.sqrt(RayTracer.trigEpsilon) <= current_t && current_t < intersecParam.intersectionDist && trigIndx != ShapeIndx){
					intersecParam.intersectionDist = current_t;
					intersecParam.closestShapeType = Shape.TRIANGLE;
					intersecParam.closestShapeIndx = trigIndx;
				}
			}
		}
		if (intersecParam.intersectionDist == Double.MAX_VALUE){
			intersecParam.closestShapeType = Shape.BACKGROUND;
			intersecParam.closestShapeIndx = -1;
		}
		 return intersecParam;

	}
	
	public static double FindLightIntersection(Ray ray,Scene scene, Shape closestShapeType, int closestShapeIndx,double tmin,Vector lightPos){
		// create a ray from camera through pixel center
				
		// check intersection with all spheres
		double totalTransparency = 1;
		if (scene.spheres != null){
			for (int sphereIndx=0; sphereIndx<scene.spheres.size();sphereIndx++){
				double current_t = sphereIntersect(ray,scene.spheres.get(sphereIndx), lightPos);
				if (0 <= current_t && current_t  < tmin ){
					Sphere blockingSphere = scene.spheres.get(sphereIndx);
					double trans = scene.materials.get((blockingSphere.sphereMatIndx)- 1).matTransparency;
					totalTransparency = totalTransparency * trans;
				}
			}
		}
		// check intersection with all planes
		if (scene.planes != null){
			for (int planeIndx=0; planeIndx<scene.planes.size();planeIndx++){
				double current_t = planeIntersect(ray,scene.planes.get(planeIndx), lightPos);
				if (0 <= current_t && current_t < tmin){
					Plane blockingPlane = scene.planes.get(planeIndx);
					double trans = scene.materials.get(blockingPlane.planeMatIndx - 1).matTransparency;
					totalTransparency = totalTransparency * trans;
				}
			}
		}
		// check intersection with all triangles
		if (scene.triangles != null){
			for (int trigIndx=0; trigIndx<scene.triangles.size();trigIndx++){
				double current_t = triangleIntersect(ray,scene.triangles.get(trigIndx),lightPos);
				if (0 <= current_t && current_t < tmin){
					Triangle blockingTrig = scene.triangles.get(trigIndx);
					double trans = scene.materials.get(blockingTrig.trigMatIndx - 1).matTransparency;
					totalTransparency = totalTransparency * trans;					
				}
			}
		}
		return totalTransparency;
	}
	
	private static double planeIntersect(Ray ray, Plane plane, Vector sourcePos) {
		if (Vector.dotProduct(ray.direction,plane.planeNorm)==0){
			return -1;
		}
		double p0dotN = Vector.dotProduct(sourcePos,plane.planeNorm);
		double VdotN = Vector.dotProduct(ray.direction,plane.planeNorm);
		double offset = plane.planeOffset;
		double t = - (p0dotN - offset) / VdotN;  
		return t;
	}
	
	////////////////////////////////////////////////////////////////////////////////////
	// Moller-Trumbore triangle intersection algorithm
	// addopted from wikipedia:
	// https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
	////////////////////////////////////////////////////////////////////////////////////
	public static double triangleIntersect(Ray ray, Triangle triangle,Vector sourcePos){
		
		Vector edge1 = Vector.subVector(triangle.trigVert2,triangle.trigVert1);
		Vector edge2 = Vector.subVector(triangle.trigVert3,triangle.trigVert1);
		Vector pVec = Vector.crossProduct(ray.direction, edge2);
		double det = Vector.dotProduct(edge1, pVec);
		if (det > -trigEpsilon && det < trigEpsilon){
			return -1;
		}
		double inversDet = 1/(double) det;
		Vector tVec = Vector.subVector(sourcePos, triangle.trigVert1);
		double u = Vector.dotProduct(tVec, pVec) * inversDet;
		if (u < 0 || u > 1){
			return -1;
		}
		Vector qVec = Vector.crossProduct(tVec, edge1);
		double v = Vector.dotProduct(ray.direction, qVec) * inversDet;
		if (v < 0 || u +v > 1){
			return -1;
		}
		double t = Vector.dotProduct(edge2, qVec) * inversDet;
		return t;
	}

	public static double sphereIntersect(Ray ray,Sphere sphere, Vector sourcePos){
		
		Vector L = Vector.subVector(sphere.spherePos, sourcePos);
		double t_ca = Vector.dotProduct(L, ray.direction);
		if (t_ca < 0){
			return -1;
		}
		double d_sqrd = Vector.dotProduct(L,L) - (t_ca*t_ca);
		double r_sqrd = (sphere.sphereRad)*(sphere.sphereRad);
		
		if (d_sqrd > r_sqrd){
			return -1;
		}
		double t_hc = Math.sqrt(r_sqrd - d_sqrd);
		double t1 = t_ca - t_hc;
		double t2 = t_ca + t_hc;
		if (t1<0){
			return t2;
		}
		return t1;
	}

	public static Vector getPixelColor(Shape closestShapeType,int closestShapeIndx , Vector IntersectPosition,Vector shapeNormal,RayTracer tracer,Random rand,Material shapeMaterial,Ray ray,int recursionDepth){
		if (recursionDepth == tracer.recMax){
			return tracer.scene.sceneSettings.backgroundRGB;
		}
		switch (closestShapeType){
			case SPHERE: shapeMaterial = tracer.scene.materials.get((tracer.scene.spheres.get(closestShapeIndx).sphereMatIndx)-1);
				 break;
			case TRIANGLE: shapeMaterial = tracer.scene.materials.get((tracer.scene.triangles.get(closestShapeIndx).trigMatIndx)-1);
				break;
			case PLANE: shapeMaterial = tracer.scene.materials.get((tracer.scene.planes.get(closestShapeIndx).planeMatIndx)-1);
				break;
			default: 
				break;
		}
		Vector Out = new Vector();
		Vector pointColor = new Vector();
		pointColor.xCor = 0; pointColor.yCor = 0; pointColor.zCor = 0;
		
		
		// sum up RGB contribution of all light sources
		// [1] loop over all lights in scene
		// [2] for each light source, check if it hits the object.
		// 	   if it does, just add the RGB light intensity
		//     if it does not, add the light intensity mult by (1-shadow intensity)
		// [3] for each light, calculate it's hit factor at the point (soft shadows)
		// 	   and mult RGB by the hitFactor

		for (int lightIndx=0; lightIndx < tracer.scene.lights.size();lightIndx++){
			// define tempPointColor
			Vector tempPointColor = new Vector();
			tempPointColor.xCor = 0; tempPointColor.yCor = 0; tempPointColor.zCor = 0;
			// check if light source hits the source directly
			Light currentLight = tracer.scene.lights.get(lightIndx);
			// create ray from light source to intersection position (with epsilon shift)
			// note that lightRay is normalized
			Vector IntersectPositionNew = Vector.addVector(IntersectPosition,Vector.scalarMult(ray.direction, -epsilon));
			Ray lightRay = new Ray(IntersectPositionNew,currentLight.lightPos);
			// calculate distance between light Source to intersection position
			// calculate soft shadows:
			// determine if there is a primitive between light source and intersection position
			//
			double hitFactor = calcHitFactor(tracer,lightRay,currentLight,rand,closestShapeType,closestShapeIndx,IntersectPosition);
			// basic color calculation
			// calculate normal
			Vector Normal = getShapeNoraml(tracer.scene,IntersectPositionNew,closestShapeType,closestShapeIndx);
			// calculate diffuse light
			Vector diffuse = diffuseLight(shapeMaterial,Vector.scalarMult(lightRay.direction,-1),Normal,tracer.scene,currentLight.lightRGB,closestShapeType,1);
			// calculate specular light
			Vector spec = specularLight(shapeMaterial,Vector.scalarMult(lightRay.direction,-1),Normal,tracer.scene,currentLight.lightRGB,1,currentLight.lightSpec,ray);
			// add spec and diffuse to pixel's color
			tempPointColor = Vector.addVector(diffuse, spec);
			tempPointColor = Vector.scalarMult(tempPointColor, hitFactor * currentLight.shadow + (1 - currentLight.shadow));
			pointColor = Vector.addVector(pointColor, tempPointColor);
		}
		//Reflection calculation
		Material mat = null;
		Vector reflaxionColor;
		Vector ref;
		Vector transparentColor;
		Vector OutRayDirection = Vector.normalizeVector(computeReflexionVector(shapeNormal,ray.direction));
		Ray OutRay = new Ray();
		OutRay.direction = OutRayDirection;
		IntersectionParameters intersecParam1 = FindIntersection(OutRay,tracer.scene,IntersectPosition,-1);
		if (Shape.BACKGROUND != intersecParam1.closestShapeType){
			    Vector IntersectPos1 = getShapePosittion(intersecParam1.intersectionDist,OutRay,IntersectPosition);
				shapeNormal = getShapeNoraml(tracer.scene,IntersectPos1,intersecParam1.closestShapeType,intersecParam1.closestShapeIndx);
				reflaxionColor = getPixelColor(intersecParam1.closestShapeType,intersecParam1.closestShapeIndx,IntersectPos1,shapeNormal,tracer,rand,mat,OutRay,recursionDepth+1);
				ref = reflaxionLight(reflaxionColor,shapeMaterial);
		}
		else {
			ref = reflaxionLight(tracer.scene.sceneSettings.backgroundRGB,shapeMaterial);	// we hit background no need to recurse
		}
		//transparency calculation
		Vector sourceEpsilon = Vector.addVector(IntersectPosition,Vector.scalarMult(ray.direction, +epsilon));
		//Vector sourceEpsilon = addEpsilon(IntersectPosition);
		IntersectionParameters intersecParam2 = FindIntersection(ray,tracer.scene,sourceEpsilon,closestShapeIndx);
		if (Shape.BACKGROUND != intersecParam2.closestShapeType){
			Vector IntesectionPos2 = getShapePosittion(intersecParam2.intersectionDist,ray,sourceEpsilon);
			shapeNormal = getShapeNoraml(tracer.scene,IntesectionPos2,intersecParam2.closestShapeType,intersecParam2.closestShapeIndx);
			mat = null;
			transparentColor = getPixelColor(intersecParam2.closestShapeType,intersecParam2.closestShapeIndx,IntesectionPos2,shapeNormal,tracer,rand,mat,ray,recursionDepth+1);
		}
		else {
			transparentColor = 	tracer.scene.sceneSettings.backgroundRGB;// we hit background no need to recurse
		}
		Vector OutTrans = Vector.scalarMult(transparentColor,shapeMaterial.matTransparency);
		Vector OutAmbDiff = Vector.scalarMult(pointColor,1 - shapeMaterial.matTransparency);
		Out = Vector.addVector(OutTrans,OutAmbDiff);
		Out = Vector.addVector(Out,ref);
		return threshhold(Out);
		
	}
	
	public static Vector addEpsilon(Vector intersection){
		return Vector.scalarMult(intersection,1.001);
	}
	
	
	public static Vector getShapeNoraml(Scene scene, Vector shapePos, Shape shapeType, int shapeIndx){
		Vector normal = new Vector();
		if (shapeType == Shape.SPHERE){
			normal = Vector.subVector(shapePos,scene.spheres.get(shapeIndx).spherePos);
			normal = Vector.normalizeVector(normal);
			return normal;
		}
		else if (shapeType == Shape.TRIANGLE){
			Vector v1 = Vector.addVector(scene.triangles.get(shapeIndx).trigVert1,Vector.scalarMult(scene.triangles.get(shapeIndx).trigVert2,-1));
			Vector v2 = Vector.addVector(scene.triangles.get(shapeIndx).trigVert2,Vector.scalarMult(scene.triangles.get(shapeIndx).trigVert3,-1));
			return Vector.normalizeVector(Vector.crossProduct(v2,v1));
		}
		else {
			return scene.planes.get(shapeIndx).planeNorm;
		}
	}
	
	static Vector threshhold(Vector v ){
		if (v.xCor >= 1){
			v.xCor = 1;
		}
		if (v.yCor >= 1){
			v.yCor = 1;
		}
		if (v.zCor >= 1){
			v.zCor = 1;
		}
		if (v.xCor < 0){
			v.xCor = 0;
		}
		if (v.yCor < 0){
			v.yCor = 0;
		}
		if (v.zCor < 0){
			v.zCor = 0;
		}
		
		return v;
	}
	
	public static Vector getShapePosittion(double intersectionDist,Ray ray,Vector Raysource){
		return Vector.addVector(Raysource, Vector.scalarMult(ray.direction, intersectionDist));
	}
	
	public static double calcHitFactor(RayTracer tracer,Ray lightRay,Light currentLight, Random rand, Shape closestShapeType,int closestShapeIndx, Vector IntersectPosition){
		// calculate light plane
		// [a] get two perpendicular vectors
		//double hitFactor = calcHitFactor(tracer,lightRay,currentLight,rand,closestShapeType,closestShapeIndx,IntersectPosition);
		Vector planeVector1 = Vector.gerArbitraryPerpendicularVector(lightRay.direction);
		if (planeVector1 == null){
			return 0;
		}
		planeVector1 = Vector.normalizeVector(planeVector1);
		Vector planeVector2 = Vector.crossProduct(lightRay.direction, planeVector1);
		planeVector2 = Vector.normalizeVector(planeVector2);
		// [b] cast many rays for each light
		int rayNum = tracer.scene.sceneSettings.rootNumSdwRays;
		double lightWidth = currentLight.lightWidth;
		double cellSize = lightWidth/rayNum;
		double rayCounter = 0;
		
		for (int Nx = 0; Nx < rayNum;Nx++){
			Vector xLocation = Vector.addVector(Vector.scalarMult(planeVector1, -lightWidth/2.0),Vector.scalarMult(planeVector1,(Nx*cellSize)));
			Vector randXshift = Vector.addVector(xLocation, Vector.scalarMult(planeVector1,cellSize*rand.nextDouble()));
			for (int Ny = 0 ; Ny <rayNum;Ny++){
				Vector yLocation = Vector.addVector(Vector.scalarMult(planeVector2, -lightWidth/2.0),Vector.scalarMult(planeVector2,(Ny*cellSize)));
				Vector randYshift = Vector.addVector(yLocation, Vector.scalarMult(planeVector2,cellSize*rand.nextDouble()));
				Vector currentLightPos = Vector.addVector(currentLight.lightPos,Vector.addVector(randXshift, randYshift));
				
				Vector IntersectPositionNew = Vector.addVector(IntersectPosition,Vector.scalarMult(lightRay.direction, -epsilon));
				Ray templightRay = new Ray(IntersectPositionNew,currentLightPos);
				// calculate distance between light Source to intersection position
				double t_dist = Math.sqrt(Vector.dotProduct(Vector.subVector(currentLightPos,IntersectPositionNew),Vector.subVector(currentLightPos,IntersectPositionNew)));
				// calculate soft shadows:
				// determine if there is a primitive between light source and intersection position
				double hitFactor = FindLightIntersection(templightRay,tracer.scene,closestShapeType,closestShapeIndx,t_dist,currentLightPos);
				rayCounter+=hitFactor;
			}
		}
		double hitFactor = rayCounter/((double) (rayNum*rayNum));
		return hitFactor;
	
	}
	
	
	public static Vector reflaxionLight(Vector reflaxionColor,Material mat){
		return Vector.elementProduct(reflaxionColor,mat.reflectionColor);
	}

	public static class RayTracerException extends Exception {
		public RayTracerException(String msg) {  super(msg); }
	}
	
	public static class Scene{
		
		Camera camera = null;
		SceneSettings sceneSettings = null;
		ArrayList<Material> materials = null;
		ArrayList<Sphere> spheres = null;
		ArrayList<Plane> planes = null;
		ArrayList<Triangle> triangles = null;
		ArrayList<Light> lights = null;
		
	}
	
	public static class Camera{
		
		Vector camPosition = new Vector();
		Vector camDirection = new Vector();
		Vector camUpVector = new Vector();
		double screenHeight;
		double screenDist;
		double screenWidth;
		
		public Camera(String[] camString){
			
			camPosition.xCor = Double.parseDouble(camString[0]);
			camPosition.yCor = Double.parseDouble(camString[1]);
			camPosition.zCor = Double.parseDouble(camString[2]);
			
			camDirection.xCor = Double.parseDouble(camString[3]);
			camDirection.yCor = Double.parseDouble(camString[4]);
			camDirection.zCor = Double.parseDouble(camString[5]);
			
			camUpVector.xCor = Double.parseDouble(camString[6]);
			camUpVector.yCor = Double.parseDouble(camString[7]);
			camUpVector.zCor = Double.parseDouble(camString[8]);
			
			screenDist = Double.parseDouble(camString[9]);
			screenWidth = Double.parseDouble(camString[10]);
			
			
		}
	}
	public static class SceneSettings{
		
		// Add code here to parse general settings parameters
		Vector backgroundRGB = new Vector();
		int rootNumSdwRays;
		int recMax;
		int samplingLevel = 1;
		
		public SceneSettings(String[] setString){
			
			backgroundRGB.xCor = Double.parseDouble(setString[0]);
			backgroundRGB.yCor = Double.parseDouble(setString[1]);
			backgroundRGB.zCor = Double.parseDouble(setString[2]);
			rootNumSdwRays = Integer.parseInt(setString[3]);
			recMax = Integer.parseInt(setString[4]);
			if (setString.length == 6){
				samplingLevel = Integer.parseInt(setString[5]);
			}
			else{
				System.out.println("No super sampling input data. Default sampling level of " + samplingLevel + " was chosen\n");
			}
		}
	}

	public static class Material {
		
		Vector difuseRGB = new Vector();
		Vector specularRGB = new Vector();
		Vector reflectionColor = new Vector();
		public double phongSpec;
		public double matTransparency;
		public double matIncidence;
		
		public Material(String[] materialString){
			
			difuseRGB.xCor = Double.parseDouble(materialString[0]);
			difuseRGB.yCor = Double.parseDouble(materialString[1]);
			difuseRGB.zCor = Double.parseDouble(materialString[2]);
			
			specularRGB.xCor = Double.parseDouble(materialString[3]);
			specularRGB.yCor = Double.parseDouble(materialString[4]);
			specularRGB.zCor = Double.parseDouble(materialString[5]);
			
			reflectionColor.xCor = Double.parseDouble(materialString[6]);
			reflectionColor.yCor = Double.parseDouble(materialString[7]);
			reflectionColor.zCor = Double.parseDouble(materialString[8]);
			
			phongSpec = Double.parseDouble(materialString[9]);
			matTransparency = Double.parseDouble(materialString[10]);
			//matIncidence = Double.parseDouble(materialString[12]);
			
		}
	}
	public static class Sphere{
		
		Vector spherePos = new Vector();
		double sphereRad;
		int sphereMatIndx; 
		
		public Sphere(String[] sphereString){
			
			spherePos.xCor = Double.parseDouble(sphereString[0]);
			spherePos.yCor = Double.parseDouble(sphereString[1]);
			spherePos.zCor = Double.parseDouble(sphereString[2]);
			
			sphereRad = Double.parseDouble(sphereString[3]);
			sphereMatIndx = Integer.parseInt(sphereString[4]);
		}
	}
	public static class Plane{
		
		Vector planeNorm = new Vector();
		double planeOffset;
		int planeMatIndx;
		
		public Plane(String[] planeString){
			
			planeNorm.xCor = Double.parseDouble(planeString[0]);
			planeNorm.yCor = Double.parseDouble(planeString[1]);
			planeNorm.zCor = Double.parseDouble(planeString[2]);
			
			planeOffset = Double.parseDouble(planeString[3]);
			planeMatIndx = Integer.parseInt(planeString[4]);
		}

		public Plane() {
		}
	}
	public static class Triangle{
		
		Vector trigVert1 = new Vector();
		Vector trigVert2 = new Vector();
		Vector trigVert3 = new Vector();
		int trigMatIndx;
		
		public Triangle(String[] triangleString){
			
			trigVert1.xCor = Double.parseDouble(triangleString[0]);
			trigVert1.yCor = Double.parseDouble(triangleString[1]);
			trigVert1.zCor = Double.parseDouble(triangleString[2]);
			
			trigVert2.xCor = Double.parseDouble(triangleString[3]);
			trigVert2.yCor = Double.parseDouble(triangleString[4]);
			trigVert2.zCor = Double.parseDouble(triangleString[5]);
			
			trigVert3.xCor = Double.parseDouble(triangleString[6]);
			trigVert3.yCor = Double.parseDouble(triangleString[7]);
			trigVert3.zCor = Double.parseDouble(triangleString[8]);
			
			trigMatIndx = Integer.parseInt(triangleString[9]);
		}
	}
	
	public static class Light{
		
		Vector lightPos = new Vector();
		Vector lightRGB = new Vector();
		
		// TODO: verify double or int
		double lightSpec;
		double shadow;
		int lightWidth;
		
		public Light(String[] lightString){
			
			lightPos.xCor = Double.parseDouble(lightString[0]);
			lightPos.yCor = Double.parseDouble(lightString[1]);
			lightPos.zCor = Double.parseDouble(lightString[2]);
			
			lightRGB.xCor = Double.parseDouble(lightString[3]);
			lightRGB.yCor = Double.parseDouble(lightString[4]);
			lightRGB.zCor = Double.parseDouble(lightString[5]);
			
			lightSpec = Double.parseDouble(lightString[6]);
			shadow = Double.parseDouble(lightString[7]);
			lightWidth = Integer.parseInt(lightString[8]);
		}
	}
	
	public static class Vector{
		
		double xCor = 0;
		double yCor = 0;
		double zCor = 0;
		
		
		public static Double dotProduct(Vector v1,Vector v2){
			return (v1.xCor*v2.xCor+v1.yCor*v2.yCor+v1.zCor*v2.zCor);
		}
		public static Vector crossProduct(Vector v1,Vector v2){
			Vector v3 = new Vector();
			v3.xCor = v1.yCor*v2.zCor - v1.zCor*v2.yCor;
			v3.yCor = -(v1.xCor*v2.zCor - v1.zCor*v2.xCor);
			v3.zCor = v1.xCor*v2.yCor - v1.yCor*v2.xCor;
			return v3;
		}
		public static Vector addVector(Vector v1,Vector v2){
			Vector v3 = new Vector();
			v3.xCor = v1.xCor + v2.xCor;
			v3.yCor = v1.yCor + v2.yCor;
			v3.zCor = v1.zCor + v2.zCor;
			return v3;
		}
		public static Vector subVector(Vector v1,Vector v2){
			Vector v3 = new Vector();
			v3.xCor = v1.xCor - v2.xCor;
			v3.yCor = v1.yCor - v2.yCor;
			v3.zCor = v1.zCor - v2.zCor;
			return v3;
		}
		public static Vector scalarMult(Vector v,double scalar){
			Vector scaledVector = new Vector();
			scaledVector.xCor = (v.xCor * scalar);
			scaledVector.yCor = (v.yCor * scalar);
			scaledVector.zCor = (v.zCor * scalar);
			return scaledVector;
		}
		public static Vector elementProduct(Vector v1, Vector v2){
			Vector v3 = new Vector();
			v3.xCor = v1.xCor * v2.xCor;
			v3.yCor = v1.yCor * v2.yCor;
			v3.zCor = v1.zCor * v2.zCor;
			return v3;
		}
		public static Vector normalizeVector(Vector v){
			
			if ((v.xCor == 0) && (v.yCor == 0) && (v.zCor == 0)){
				return v;
			}
			Vector vNorm = new Vector();
			double vSize = Math.sqrt(Vector.dotProduct(v,v));
			vNorm = Vector.scalarMult(v,(1/vSize));
			return vNorm;
		}
		
		// TODO: if null is returned, it means that the inpur was the zero vector
		// need to treat this situation
		public static Vector gerArbitraryPerpendicularVector(Vector v){
			
			Vector perpVector = new Vector();
			if (v.xCor == 0 && v.yCor == 0){
				if (v.zCor == 0){
					return null;
				}
				perpVector.xCor = 0;
				perpVector.yCor = 1;
				perpVector.zCor = 0;
				return perpVector;
			}
			perpVector.xCor = -v.yCor;
			perpVector.yCor = v.xCor;
			perpVector.zCor = 0;
			return perpVector;
		}
		
		public static Vector decreaseEpsilon(Vector v){
			return Vector.scalarMult(v, 1-RayTracer.epsilon);
			
		}
		
	}	

	public static class Ray{
		Vector direction;
		
		public Ray() {}
		Ray(Vector source, Vector target){
			this.direction = Vector.subVector(source,target);
			this.direction = Vector.normalizeVector(this.direction);
		}
		
		Ray(RayTracer tracer, Scene scene,int i,int j,int Ni,int Nj,Random rand,int imageHeight,int imageWidth){
			
			scene.camera.screenHeight = scene.camera.screenWidth*imageHeight/imageWidth;
			Vector towards = Vector.normalizeVector(Vector.addVector(scene.camera.camDirection,Vector.scalarMult(scene.camera.camPosition,-1)));
			Vector left = Vector.normalizeVector(Vector.crossProduct(scene.camera.camUpVector,towards));
			Vector up = Vector.normalizeVector(Vector.crossProduct(left,towards));
			Vector lowerLeftCorner = Vector.addVector(scene.camera.camPosition,
					Vector.addVector(Vector.scalarMult(towards,scene.camera.screenDist)
					,Vector.addVector(Vector.scalarMult(left,-0.5*scene.camera.screenWidth)
					,Vector.scalarMult(up,-0.5*scene.camera.screenHeight))));
			Vector currentP = Vector.addVector(lowerLeftCorner,Vector.
					addVector(Vector.scalarMult(left,i*(scene.camera.screenWidth/imageHeight)),
					Vector.scalarMult(up,j*(scene.camera.screenWidth/imageWidth))));
			
			double pixelWidth = scene.camera.screenWidth/ ((double) tracer.imageWidth);
			double pixelHeight = scene.camera.screenHeight/ ((double) tracer.imageHeight);
			
			double randomLeft, randomUp;
			do {
				randomLeft = rand.nextDouble();
			} while(randomLeft == 0);
			randomLeft = randomLeft - 0.5; 
			do {
				randomUp = rand.nextDouble();
			} while(randomUp == 0);
			randomUp = randomUp -0.5;
			int samp = scene.sceneSettings.samplingLevel;
			if (samp == 1){
				randomLeft = 0; randomUp=0;
			}
			Vector currentShift = Vector.addVector(
					Vector.scalarMult(left,  (Ni+randomLeft)*pixelWidth/(double)samp), 
					Vector.scalarMult(up,  (Nj+randomUp)*pixelHeight/(double)samp));
			currentP = Vector.addVector(currentP,currentShift);
			direction = Vector.normalizeVector(Vector.subVector(currentP, scene.camera.camPosition));
			
		}
	}
	
	public static class IntersectionParameters{
		
		double intersectionDist; 
		Shape closestShapeType; 
		int closestShapeIndx;
			
	}

	//Vector diffuse = diffuseLight(shapeMaterial,lightRay.direction,Normal,tracer.scene,currentLight.lightRGB,1);
	public static Vector diffuseLight(Material mat,Vector L, Vector Normal,Scene scene,Vector Ip,Shape ShapeType,double hitFactor){
		Vector Out = new Vector();
		double prod;
		if (ShapeType == Shape.SPHERE || ShapeType == Shape.PLANE){
			prod = (Vector.dotProduct(Vector.normalizeVector(L),Vector.normalizeVector(Normal)));
		}
		else{
			prod = java.lang.Math.abs(Vector.dotProduct(Vector.normalizeVector(L),Vector.normalizeVector(Normal)));
		}
		if (prod < 0){
			return Out;
		}
		Out = Vector.scalarMult(Vector.elementProduct(mat.difuseRGB, Ip),prod*hitFactor);
		return Out;
	}
	
	public static Vector specularLight(Material mat,Vector L, Vector Normal,Scene scene,Vector Ip,double hitFactor,double LightSpec,Ray ray){
		Vector Out = new Vector();
		Vector reflaxionVector = computeReflexionVector(Normal,ray.direction);
		double prod = java.lang.Math.abs(Vector.dotProduct(reflaxionVector,L));
		prod = Vector.dotProduct(reflaxionVector,L);
		if (prod < 0){
		return Out;
		}
		Out = Vector.scalarMult(Vector.elementProduct(mat.specularRGB,Ip),Math.pow(prod,mat.phongSpec)*hitFactor*LightSpec);
		return threshhold(Out);
	}

	public static Vector computeReflexionVector(Vector Normal,Vector v){
		double prod = -2*Vector.dotProduct(v,Normal)/Vector.dotProduct(Normal,Normal);
		Vector n = Vector.scalarMult(Normal,prod);
		return Vector.addVector(v,n); 
	}
}

	