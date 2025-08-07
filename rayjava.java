import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

// Basic 3D Vector class
class Vec3f {
    public float x, y, z;

    public Vec3f(float x, float y, float z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public Vec3f() {
        this(0, 0, 0);
    }

    // Dot product
    public float dot(Vec3f v) {
        return this.x * v.x + this.y * v.y + this.z * v.z;
    }

    // Scalar multiplication
    public Vec3f multiply(float s) {
        return new Vec3f(this.x * s, this.y * s, this.z * s);
    }

    // Vector addition
    public Vec3f add(Vec3f v) {
        return new Vec3f(this.x + v.x, this.y + v.y, this.z + v.z);
    }

    // Vector subtraction
    public Vec3f subtract(Vec3f v) {
        return new Vec3f(this.x - v.x, this.y - v.y, this.z - v.z);
    }

    // Element-wise multiplication (for color blending)
    public Vec3f hadamard(Vec3f v) {
        return new Vec3f(this.x * v.x, this.y * v.y, this.z * v.z);
    }

    // Magnitude (norm)
    public float norm() {
        return (float) Math.sqrt(this.dot(this));
    }

    // Normalize (return a new normalized vector)
    public Vec3f normalize() {
        float n = norm();
        return n > 0 ? new Vec3f(x / n, y / n, z / n) : new Vec3f(0,0,0);
    }

    // Clamping for color output
    public float get(int index) {
        switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: throw new IndexOutOfBoundsException("Invalid index for Vec3f");
        }
    }
}

// Basic 4D Vector class (for albedo)
class Vec4f {
    public float x, y, z, w;

    public Vec4f(float x, float y, float z, float w) {
        this.x = x;
        this.y = y;
        this.z = z;
        this.w = w;
    }

    public Vec4f(float val) {
        this(val, val, val, val);
    }

    public Vec4f() {
        this(0, 0, 0, 0);
    }

    public float get(int index) {
        switch (index) {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            case 3: return w;
            default: throw new IndexOutOfBoundsException("Invalid index for Vec4f");
        }
    }
}

class Light {
    public Vec3f position;
    public float intensity;

    public Light(Vec3f p, float i) {
        this.position = p;
        this.intensity = i;
    }
}

class Material {
    public float refractiveIndex;
    public Vec4f albedo;
    public Vec3f diffuseColor;
    public float specularExponent;

    public Material(float r, Vec4f a, Vec3f color, float spec) {
        this.refractiveIndex = r;
        this.albedo = a;
        this.diffuseColor = color;
        this.specularExponent = spec;
    }

    public Material() {
        this(1, new Vec4f(1, 0, 0, 0), new Vec3f(), 0);
    }
}

class Sphere {
    public Vec3f center;
    public float radius;
    public Material material;

    public Sphere(Vec3f c, float r, Material m) {
        this.center = c;
        this.radius = r;
        this.material = m;
    }

    // Returns t0, or Float.MAX_VALUE if no intersection
    public float rayIntersect(Vec3f orig, Vec3f dir) {
        Vec3f L = center.subtract(orig);
        float tca = L.dot(dir);
        float d2 = L.dot(L) - tca * tca;
        if (d2 > radius * radius) return Float.MAX_VALUE;
        float thc = (float) Math.sqrt(radius * radius - d2);
        float t0 = tca - thc;
        float t1 = tca + thc;

        if (t0 < 0) t0 = t1;
        if (t0 < 0) return Float.MAX_VALUE;
        return t0;
    }
}

public class RayTracer {

    private static Vec3f reflect(Vec3f I, Vec3f N) {
        return I.subtract(N.multiply(2.0f * I.dot(N)));
    }

    // Snell's law
    private static Vec3f refract(Vec3f I, Vec3f N, float refractiveIndex) {
        float cosi = -Math.max(-1.0f, Math.min(1.0f, I.dot(N)));
        float etai = 1, etat = refractiveIndex;
        Vec3f n = N;
        if (cosi < 0) { // if the ray is inside the object, swap the indices and invert the normal to get the correct result
            cosi = -cosi;
            float temp = etai;
            etai = etat;
            etat = temp;
            n = N.multiply(-1.0f);
        }
        float eta = etai / etat;
        float k = 1 - eta * eta * (1 - cosi * cosi);
        return k < 0 ? new Vec3f(0, 0, 0) : I.multiply(eta).add(n.multiply(eta * cosi - (float) Math.sqrt(k)));
    }

    private static boolean sceneIntersect(Vec3f orig, Vec3f dir, List<Sphere> spheres,
                                          Vec3f[] hit, Vec3f[] N, Material[] material) {
        float spheresDist = Float.MAX_VALUE;
        int hitSphereIdx = -1;

        for (int i = 0; i < spheres.size(); i++) {
            float distI = spheres.get(i).rayIntersect(orig, dir);
            if (distI != Float.MAX_VALUE && distI < spheresDist) {
                spheresDist = distI;
                hitSphereIdx = i;
            }
        }

        float checkerboardDist = Float.MAX_VALUE;
        if (Math.abs(dir.y) > 1e-3) {
            float d = -(orig.y + 4) / dir.y; // the checkerboard plane has equation y = -4
            Vec3f pt = orig.add(dir.multiply(d));
            if (d > 0 && Math.abs(pt.x) < 10 && pt.z < -10 && pt.z > -30 && d < spheresDist) {
                checkerboardDist = d;
                hit[0] = pt;
                N[0] = new Vec3f(0, 1, 0);
                // Create a temporary material for the checkerboard
                material[0] = new Material(1.0f, new Vec4f(1,0,0,0), new Vec3f(), 0); // Placeholder albedo and spec
                if (((int)(0.5 * hit[0].x + 1000) + (int)(0.5 * hit[0].z)) % 2 == 1) {
                    material[0].diffuseColor = new Vec3f(0.3f, 0.3f, 0.3f);
                } else {
                    material[0].diffuseColor = new Vec3f(0.3f, 0.2f, 0.1f);
                }
                 // Set albedo components for checkerboard (diffuse only)
                material[0].albedo = new Vec4f(1.0f, 0.0f, 0.0f, 0.0f); // Pure diffuse for checkerboard
                material[0].specularExponent = 0.0f; // No specular for checkerboard
            }
        }

        if (hitSphereIdx != -1 && spheresDist < checkerboardDist) {
            hit[0] = orig.add(dir.multiply(spheresDist));
            N[0] = (hit[0].subtract(spheres.get(hitSphereIdx).center)).normalize();
            material[0] = spheres.get(hitSphereIdx).material;
        }

        return Math.min(spheresDist, checkerboardDist) < 1000;
    }


    private static Vec3f castRay(Vec3f orig, Vec3f dir, List<Sphere> spheres, List<Light> lights, int depth) {
        // Using arrays for mutable Vec3f and Material references
        Vec3f[] point = {new Vec3f()};
        Vec3f[] N = {new Vec3f()};
        Material[] material = {new Material()};

        if (depth > 4 || !sceneIntersect(orig, dir, spheres, point, N, material)) {
            return new Vec3f(0.2f, 0.7f, 0.8f); // background color
        }

        Vec3f reflectDir = reflect(dir, N[0]).normalize();
        Vec3f refractDir = refract(dir, N[0], material[0].refractiveIndex).normalize();

        // offset the original point to avoid occlusion by the object itself
        Vec3f reflectOrig = reflectDir.dot(N[0]) < 0 ? point[0].subtract(N[0].multiply(1e-3f)) : point[0].add(N[0].multiply(1e-3f));
        Vec3f refractOrig = refractDir.dot(N[0]) < 0 ? point[0].subtract(N[0].multiply(1e-3f)) : point[0].add(N[0].multiply(1e-3f));

        Vec3f reflectColor = castRay(reflectOrig, reflectDir, spheres, lights, depth + 1);
        Vec3f refractColor = castRay(refractOrig, refractDir, spheres, lights, depth + 1);

        float diffuseLightIntensity = 0;
        float specularLightIntensity = 0;

        for (Light light : lights) {
            Vec3f lightDir = (light.position.subtract(point[0])).normalize();
            float lightDistance = (light.position.subtract(point[0])).norm();

            // checking if the point lies in the shadow of the light
            Vec3f shadowOrig = lightDir.dot(N[0]) < 0 ? point[0].subtract(N[0].multiply(1e-3f)) : point[0].add(N[0].multiply(1e-3f));
            Vec3f[] shadowPt = {new Vec3f()};
            Vec3f[] shadowN = {new Vec3f()};
            Material[] tmpMaterial = {new Material()};

            if (sceneIntersect(shadowOrig, lightDir, spheres, shadowPt, shadowN, tmpMaterial) &&
                (shadowPt[0].subtract(shadowOrig)).norm() < lightDistance) {
                continue; // In shadow, skip this light
            }

            diffuseLightIntensity += light.intensity * Math.max(0.0f, lightDir.dot(N[0]));
            specularLightIntensity += Math.pow(Math.max(0.0f, -reflect(lightDir.multiply(-1.0f), N[0]).dot(dir)), material[0].specularExponent) * light.intensity;
        }

        // Combine all light components and material properties
        Vec3f finalColor = material[0].diffuseColor.hadamard(new Vec3f(diffuseLightIntensity, diffuseLightIntensity, diffuseLightIntensity)).multiply(material[0].albedo.get(0));
        finalColor = finalColor.add(new Vec3f(1.0f, 1.0f, 1.0f).multiply(specularLightIntensity * material[0].albedo.get(1)));
        finalColor = finalColor.add(reflectColor.multiply(material[0].albedo.get(2)));
        finalColor = finalColor.add(refractColor.multiply(material[0].albedo.get(3)));

        return finalColor;
    }

    private static void render(List<Sphere> spheres, List<Light> lights) {
        final int width = 1024;
        final int height = 768;
        final float fov = (float) (Math.PI / 2.0);
        Vec3f[] framebuffer = new Vec3f[width * height];

        // Java does not have a direct equivalent for OpenMP's #pragma omp parallel for
        // For actual parallelism, you would use ExecutorService or parallel streams.
        // For simplicity, this is a sequential loop.
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                float x = (2 * (i + 0.5f) / (float) width - 1) * (float) Math.tan(fov / 2.0f) * width / (float) height;
                float y = -(2 * (j + 0.5f) / (float) height - 1) * (float) Math.tan(fov / 2.0f);
                Vec3f dir = new Vec3f(x, y, -1).normalize();
                framebuffer[i + j * width] = castRay(new Vec3f(0, 0, 0), dir, spheres, lights, 0);
            }
        }

        String filename = "./out.ppm";
        try (FileOutputStream ofs = new FileOutputStream(filename)) {
            ofs.write(("P6\n" + width + " " + height + "\n255\n").getBytes());

            for (int i = 0; i < height * width; ++i) {
                Vec3f c = framebuffer[i];
                float max = Math.max(c.x, Math.max(c.y, c.z));
                if (max > 1) {
                    c = c.multiply(1.0f / max);
                }

                ofs.write((byte) (255 * Math.max(0.0f, Math.min(1.0f, c.x))));
                ofs.write((byte) (255 * Math.max(0.0f, Math.min(1.0f, c.y))));
                ofs.write((byte) (255 * Math.max(0.0f, Math.min(1.0f, c.z))));
            }
        } catch (IOException e) {
            System.err.println("Error writing PPM file: " + e.getMessage());
        }
        System.out.println("Rendered image to " + filename);
    }

    public static void main(String[] args) {
        Material ivory = new Material(1.0f, new Vec4f(0.6f, 0.3f, 0.1f, 0.0f), new Vec3f(0.4f, 0.4f, 0.3f), 50.0f);
        Material glass = new Material(1.5f, new Vec4f(0.0f, 0.5f, 0.1f, 0.8f), new Vec3f(0.6f, 0.7f, 0.8f), 125.0f);
        Material redRubber = new Material(1.0f, new Vec4f(0.9f, 0.1f, 0.0f, 0.0f), new Vec3f(0.3f, 0.1f, 0.1f), 10.0f);
        Material mirror = new Material(1.0f, new Vec4f(0.0f, 10.0f, 0.8f, 0.0f), new Vec3f(1.0f, 1.0f, 1.0f), 1425.0f);

        List<Sphere> spheres = new ArrayList<>();
        spheres.add(new Sphere(new Vec3f(-3, 0, -16), 2, glass));
        spheres.add(new Sphere(new Vec3f(-1.0f, -1.5f, -12), 2, glass));
        spheres.add(new Sphere(new Vec3f(1.5f, -0.5f, -18), 3, redRubber));
        spheres.add(new Sphere(new Vec3f(7, 5, -18), 4, mirror));

        List<Light> lights = new ArrayList<>();
        lights.add(new Light(new Vec3f(-20, 20, 20), 1.5f));
        lights.add(new Light(new Vec3f(30, 50, -25), 1.8f));
        lights.add(new Light(new Vec3f(30, 20, 30), 1.7f));

        render(spheres, lights);
    }
}