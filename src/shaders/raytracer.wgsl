const THREAD_COUNT = 16;
const RAY_TMIN = 0.0001;
const RAY_TMAX = 100.0;
const PI = 3.1415927f;
const FRAC_1_PI = 0.31830987f;
const FRAC_2_PI = 1.5707964f;

@group(0) @binding(0)  
  var<storage, read_write> fb : array<vec4f>;

@group(0) @binding(1)
  var<storage, read_write> rtfb : array<vec4f>;

@group(1) @binding(0)
  var<storage, read_write> uniforms : array<f32>;

@group(2) @binding(0)
  var<storage, read_write> spheresb : array<sphere>;

@group(2) @binding(1)
  var<storage, read_write> quadsb : array<quad>;

@group(2) @binding(2)
  var<storage, read_write> boxesb : array<box>;

@group(2) @binding(3)
  var<storage, read_write> trianglesb : array<triangle>;

@group(2) @binding(4)
  var<storage, read_write> meshb : array<mesh>;

struct ray {
  origin : vec3f,
  direction : vec3f,
};

struct sphere {
  transform : vec4f,
  color : vec4f,
  material : vec4f,
};

struct quad {
  Q : vec4f,
  u : vec4f,
  v : vec4f,
  color : vec4f,
  material : vec4f,
};

struct box {
  center : vec4f,
  radius : vec4f,
  rotation: vec4f,
  color : vec4f,
  material : vec4f,
};

struct triangle {
  v0 : vec4f,
  v1 : vec4f,
  v2 : vec4f,
};

struct mesh {
  transform : vec4f,
  scale : vec4f,
  rotation : vec4f,
  color : vec4f,
  material : vec4f,
  min : vec4f,
  max : vec4f,
  show_bb : f32,
  start : f32,
  end : f32,
};

struct material_behaviour {
  scatter : bool,
  direction : vec3f,
};

struct camera {
  origin : vec3f,
  lower_left_corner : vec3f,
  horizontal : vec3f,
  vertical : vec3f,
  u : vec3f,
  v : vec3f,
  w : vec3f,
  lens_radius : f32,
};

struct hit_record {
  t : f32,
  p : vec3f,
  normal : vec3f,
  object_color : vec4f,
  object_material : vec4f,
  frontface : bool,
  hit_anything : bool,
};

fn ray_at(r: ray, t: f32) -> vec3f
{
  return r.origin + t * r.direction;
}

fn get_ray(cam: camera, uv: vec2f, rng_state: ptr<function, u32>) -> ray
{
  var rd = cam.lens_radius * rng_next_vec3_in_unit_disk(rng_state);
  var offset = cam.u * rd.x + cam.v * rd.y;
  return ray(cam.origin + offset, normalize(cam.lower_left_corner + uv.x * cam.horizontal + uv.y * cam.vertical - cam.origin - offset));
}

fn get_camera(lookfrom: vec3f, lookat: vec3f, vup: vec3f, vfov: f32, aspect_ratio: f32, aperture: f32, focus_dist: f32) -> camera
{
  var camera = camera();
  camera.lens_radius = aperture / 2.0;

  var theta = degrees_to_radians(vfov);
  var h = tan(theta / 2.0);
  var w = aspect_ratio * h;

  camera.origin = lookfrom;
  camera.w = normalize(lookfrom - lookat);
  camera.u = normalize(cross(vup, camera.w));
  camera.v = cross(camera.u, camera.w);

  camera.lower_left_corner = camera.origin - w * focus_dist * camera.u - h * focus_dist * camera.v - focus_dist * camera.w;
  camera.horizontal = 2.0 * w * focus_dist * camera.u;
  camera.vertical = 2.0 * h * focus_dist * camera.v;

  return camera;
}

fn environment_color(direction: vec3f, color1: vec3f, color2: vec3f) -> vec3f
{
  var unit_direction = normalize(direction);
  var t = 0.5 * (unit_direction.y + 1.0);
  var col = (1.0 - t) * color1 + t * color2;

  var sun_direction = normalize(vec3(uniforms[13], uniforms[14], uniforms[15]));
  var sun_color = int_to_rgb(i32(uniforms[17]));
  var sun_intensity = uniforms[16];
  var sun_size = uniforms[18];

  var sun = clamp(dot(sun_direction, unit_direction), 0.0, 1.0);
  col += sun_color * max(0, (pow(sun, sun_size) * sun_intensity));

  return col;
}

fn check_ray_collision(r: ray, max: f32) -> hit_record
{
  var spheresCount = i32(uniforms[19]);
  var quadsCount = i32(uniforms[20]);
  var boxesCount = i32(uniforms[21]);
  var trianglesCount = i32(uniforms[22]);
  var meshCount = i32(uniforms[27]);

  var record = hit_record(max, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);

  for(var i = 0; i < spheresCount; i++) {
    var record_sphere = hit_record(max, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    var sphere_i = spheresb[i];
    var sphere_center = vec3f(sphere_i.transform[0],sphere_i.transform[1], sphere_i.transform[2]);
    var sphere_radius = f32(sphere_i.transform[3]);
    hit_sphere(sphere_center, sphere_radius, r, &record_sphere, max);

    if (record_sphere.hit_anything && record_sphere.t < record.t) {
      record = record_sphere;
      record.object_color = sphere_i.color;
      record.object_material = sphere_i.material;
    }
  }

  for (var i = 0; i < boxesCount; i++) {
    var record_box = hit_record(max, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    var box_i = boxesb[i];
    var box_center = vec3f(box_i.center[0], box_i.center[1], box_i.center[2]);
    var box_radius = vec3f(box_i.radius[0], box_i.radius[1], box_i.radius[2]);

    var q = quaternion_from_euler(vec3f(box_i.rotation[0], box_i.rotation[1], box_i.rotation[2]));
    var q_inv = q_inverse(q);

    var r_rotated = rotate_ray_quaternion(r, box_center, q);

    if (box_radius[0] < 0) {
      hit_column(r_rotated, box_center, box_radius, &record_box, max);
    }
    else {
       hit_box(r_rotated, box_center, box_radius, &record_box, max);
    }

    if (record_box.hit_anything && record_box.t < record.t) {
      record = record_box;
      var p_vec = record.p - box_center;
      record.p = box_center + rotate_vector(p_vec, q_inv);
      record.normal = rotate_vector(record.normal, q_inv);
      record.object_color = box_i.color;
      record.object_material = box_i.material;
    }

  }

  for (var i = 0; i < quadsCount; i++) {
    var record_quad = hit_record(max, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
    var quad_i = quadsb[i];
    hit_quad(r, quad_i.Q, quad_i.u, quad_i.v, &record_quad, max);

    if (record_quad.hit_anything && record_quad.t < record.t) {
      record = record_quad;
      record.object_color = quad_i.color;
      record.object_material = quad_i.material;
    }
  }

  for (var i = 0; i < meshCount; i++) {
    var mesh_i = meshb[i];
    for (var t = mesh_i.start; t <= mesh_i.end; t += 1.0) {
      var record_triangle = hit_record(max, vec3f(0.0), vec3f(0.0), vec4f(0.0), vec4f(0.0), false, false);
      var triangle_t = trianglesb[u32(t)];
      var v0 = vec4f(triangle_t.v0[0], triangle_t.v0[1], triangle_t.v0[2], 1.0);
      var v1 = vec4f(triangle_t.v1[0], triangle_t.v1[1], triangle_t.v1[2], 1.0);
      var v2 = vec4f(triangle_t.v2[0], triangle_t.v2[1], triangle_t.v2[2], 1.0);

      // scale:
      var matriz_scale: mat4x4<f32> = mat4x4<f32>(
      vec4f(mesh_i.scale[0], 0.0, 0.0, 0.0), // Coluna 1
      vec4f(0.0, mesh_i.scale[1], 0.0, 0.0), // Coluna 2
      vec4f(0.0, 0.0, mesh_i.scale[2], 0.0), // Coluna 3
      vec4f(0.0, 0.0, 0.0, 1.0)  );

      // transform:
      var matriz_transform: mat4x4<f32> = mat4x4<f32>(
      vec4f(1.0, 0.0, 0.0, 0.0), // Coluna 1
      vec4f(0.0, 1.0, 0.0, 0.0), // Coluna 2
      vec4f(0.0, 0.0, 1.0, 0.0), // Coluna 3
      vec4f(mesh_i.transform[0], mesh_i.transform[1], mesh_i.transform[2], 1.0)  );

      // rotate:
      var q = quaternion_from_euler(vec3f(-mesh_i.rotation[0], -mesh_i.rotation[1], -mesh_i.rotation[2]));
      var matriz_rotate: mat4x4<f32> = mat4x4<f32>(
      vec4f(1.0 - 2*(pow(q.y, 2)+pow(q.z, 2)), 2*(q.x*q.y + q.z*q.w), 2*(q.x*q.z - q.y*q.w), 0.0), // Coluna 1
      vec4f(2.0*(q.x*q.y - q.z*q.w), 1.0 - 2*(pow(q.x, 2)+pow(q.z, 2)), 2*(q.y*q.z + q.x*q.w), 0.0), // Coluna 2
      vec4f(2.0*(q.x*q.z + q.y*q.w), 2*(q.y*q.z - q.x*q.w), 1.0 - 2*(pow(q.x, 2)+pow(q.y, 2)), 0.0), // Coluna 3
      vec4f(0.0, 0.0, 0.0, 1.0)  );

      v0 = matriz_transform * matriz_rotate * matriz_scale * v0;
      v1 = matriz_transform * matriz_rotate * matriz_scale * v1;
      v2 = matriz_transform * matriz_rotate * matriz_scale * v2;

      var v0_t = vec3f(v0[0], v0[1], v0[2]);
      var v1_t = vec3f(v1[0], v1[1], v1[2]);
      var v2_t = vec3f(v2[0], v2[1], v2[2]);

      hit_triangle(r, v0_t, v1_t, v2_t, &record_triangle, max);

      if (record_triangle.hit_anything && record_triangle.t < record.t) {
        record = record_triangle;
        record.object_color = mesh_i.color;
        record.object_material = mesh_i.material;
      }
    }
  }

  var closest = record;
  return closest;
}

fn lambertian(normal : vec3f, absorption: f32, random_sphere: vec3f, rng_state: ptr<function, u32>) -> material_behaviour
{
  var random_value = rng_next_float(rng_state);
  if (random_value < absorption) {
    return material_behaviour(false, vec3f(0.0));
  }
  else {
    var new_direction = normalize(normal + random_sphere);
    return material_behaviour(true, new_direction);
  }
}

fn metal(normal : vec3f, direction: vec3f, fuzz: f32, random_sphere: vec3f) -> material_behaviour
{
  var new_direction = normalize(reflect(direction, normal) + fuzz*random_sphere);
  return material_behaviour(true, new_direction);
}

fn dielectric(normal : vec3f, r_direction: vec3f, refraction_index: f32, frontface: bool, random_sphere: vec3f, fuzz: f32, rng_state: ptr<function, u32>) -> material_behaviour
{ 
  var n1 = f32(1.0);
  var n2 = refraction_index;
  var n = -normal;
  if (!frontface) {
    n1 = refraction_index;
    n2 = f32(1.0);
    n = normal;
  }

  var cosseno = dot(r_direction, n);
  var seno = sqrt(1-pow(cosseno,2));
  var r0 = pow((n1-n2)/(n1+n2), 2);
  var schlick = r0 + (1-r0)*pow((1-cosseno),5);
  var new_direction = vec3f(0.0);
  var random_value = rng_next_float(rng_state);

  if (n1*seno/n2 > 1 || random_value < schlick) { // reflexão
    new_direction = normalize(reflect(r_direction, n));
  }
  else { // refração
    var perpendicular = (n1/n2)*(r_direction + dot(-r_direction, normal)*normal);
    var modulo2 = pow(perpendicular[0],2) + pow(perpendicular[1],2) + pow(perpendicular[2],2);
    var paralelo = -sqrt(1-modulo2)*normal;
    new_direction = normalize(paralelo + perpendicular);
  }

  return material_behaviour(true, new_direction);
}

fn emmisive(color: vec3f, light: f32) -> material_behaviour
{
  return material_behaviour(false, vec3f(0.0));
}

fn trace(r: ray, rng_state: ptr<function, u32>) -> vec3f
{
  var maxbounces = i32(uniforms[2]);
  var light = vec3f(0.0);
  var color = vec3f(1.0);
  var r_ = r;

  var backgroundcolor1 = int_to_rgb(i32(uniforms[11]));
  var backgroundcolor2 = int_to_rgb(i32(uniforms[12]));

  for (var j = 0; j < maxbounces; j = j + 1)
  {
    // pega objeto mais próximo onde houve colisão
    var closest = check_ray_collision(r_, RAY_TMAX);

    if (closest.hit_anything) {
      var object_color = vec3(closest.object_color[0], closest.object_color[1], closest.object_color[2]);
      var object_material = closest.object_material;

      // propriedades do material
      var smoothness = object_material[0];
      var absortion = object_material[1];
      var fuzz = object_material[1];
      var specular = object_material[2];
      var refraction_index = object_material[2];
      var light_object = object_material[3];

      var random_direction = rng_next_vec3_in_unit_sphere(rng_state);

      // Ajustando color de acordo com o material
      if (smoothness < 0) { // material dielétrico
        var dielectric_behavior = dielectric(closest.normal, r_.direction, refraction_index, closest.frontface, random_direction, fuzz, rng_state);
        r_.direction = dielectric_behavior.direction;
        color *= object_color;
      }
      else {
        var lambertian_behavior = lambertian(closest.normal, absortion, random_direction, rng_state);
        var metal_behavior = metal(closest.normal, r_.direction, fuzz, random_direction);
        if (specular == 0) {
          if (smoothness == 0) { // material lambertiano
            color *= object_color;
            if (lambertian_behavior.scatter) {
              r_.direction = lambertian_behavior.direction;
            }
            else { // raio foi absorvido
              light += environment_color(r_.direction, backgroundcolor1, backgroundcolor2) * color;
              return light;
            }
          }
          else if (smoothness > 0) { // material metálico
            r_.direction = metal_behavior.direction;
            color *= smoothness*vec3(1.0) + (1 - smoothness)*object_color;
          }
        }
        else { // material especular
          var random_float = rng_next_float(rng_state);
          if (random_float > specular) {
            color *= object_color;
            if (lambertian_behavior.scatter) {
              r_.direction = lambertian_behavior.direction;
            }
            else { // raio foi absorvido
              light += environment_color(r_.direction, backgroundcolor1, backgroundcolor2) * color;
              return light;
            }
          }
          else {
            r_.direction = metal_behavior.direction;
            color *= smoothness*vec3(1.0) + (1 - smoothness)*object_color;
          }
        }
      }

      // ajustando light em materiais emissivos
      if (light_object > 0) {
        var adicional_light = light_object*object_color;
        light += adicional_light*color;
      }
      
      r_.origin = closest.p;
    }
    else {
      light += environment_color(r_.direction, backgroundcolor1, backgroundcolor2) * color;
      return light;
    }
  }
  light += environment_color(r_.direction, backgroundcolor1, backgroundcolor2) * color;
  return light;
}

@compute @workgroup_size(THREAD_COUNT, THREAD_COUNT, 1)
fn render(@builtin(global_invocation_id) id : vec3u)
{
    var rez = uniforms[1];
    var time = u32(uniforms[0]);

    // init_rng (random number generator) we pass the pixel position, resolution and frame
    var rng_state = init_rng(vec2(id.x, id.y), vec2(u32(rez)), time);

    // Get uv
    var fragCoord = vec2f(f32(id.x), f32(id.y));
    var uv = (fragCoord + sample_square(&rng_state)) / vec2(rez);

    // Camera
    var lookfrom = vec3(uniforms[7], uniforms[8], uniforms[9]);
    var lookat = vec3(uniforms[23], uniforms[24], uniforms[25]);

    // Get camera
    var cam = get_camera(lookfrom, lookat, vec3(0.0, 1.0, 0.0), uniforms[10], 1.0, uniforms[6], uniforms[5]);
    var samples_per_pixel = i32(uniforms[4]);

    //var color = vec3(rng_next_float(&rng_state), rng_next_float(&rng_state), rng_next_float(&rng_state));
    var color = vec3f(0.0);

    // Steps:
    // 1. Loop for each sample per pixel
    // 2. Get ray
    // 3. Call trace function
    // 4. Average the color

    for (var i = 0; i < samples_per_pixel; i++) {
      var ray = get_ray(cam, uv, &rng_state);
      var light = trace(ray, &rng_state);
      color += light;
    }

    color = color/f32(samples_per_pixel);

    var color_out = vec4(linear_to_gamma(color), 1.0);
    var map_fb = mapfb(id.xy, rez);

    // 5. Accumulate the color
    var should_accumulate = uniforms[3];
    var accumulated_color = rtfb[map_fb] * should_accumulate + color_out;

    // Set the color to the framebuffer
    rtfb[map_fb] = accumulated_color;
    fb[map_fb] = accumulated_color / rtfb[map_fb].w;

    //var should_accumulate = uniforms[3];

    //rtfb[map_fb] = color_out;
    //fb[map_fb] = color_out;
}