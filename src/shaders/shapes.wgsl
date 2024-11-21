fn hit_sphere(center: vec3f, radius: f32, r: ray, record: ptr<function, hit_record>, max: f32)
{
  var a = dot(r.direction, r.direction);
  var b = 2*dot(r.direction, (r.origin - center));
  var c = dot((r.origin - center), (r.origin - center)) - radius * radius;;
  var delta = pow(b,2) - 4*a*c;

  var t = f32(0.0);
  if (delta < 0) {
    record.hit_anything = false;
    return;
  }
  if (delta > 0) {
    var t1 = (-b + sqrt(delta))/(2*a);
    var t2 = (-b - sqrt(delta))/(2*a);
    if (t1 < t2 && t1 > RAY_TMIN) {
      t = t1;
    }
    else {
      t = t2;
    }
  }
  if (delta == 0) {
    t = -b/(2*a);
  }
  if (t < RAY_TMIN || t > max) {
      record.hit_anything = false;
      return;
  }

  var intersection = ray_at(r,t);
  var n = normalize(intersection - center);
  var normal = n;

  record.t = t;
  record.p = intersection;
  record.normal = normal;
  record.hit_anything = true;

  if (dot(r.direction, normal) > 0) {
    record.frontface = false;
  }
  else {
    record.frontface = true;
  }
}

fn hit_quad(r: ray, Q: vec4f, u: vec4f, v: vec4f, record: ptr<function, hit_record>, max: f32)
{
  var n = cross(u.xyz, v.xyz);
  var normal = normalize(n);
  var D = dot(normal, Q.xyz);
  var w = n / dot(n.xyz, n.xyz);

  var denom = dot(normal, r.direction);
  if (abs(denom) < 0.0001)
  {
    record.hit_anything = false;
    return;
  }

  var t = (D - dot(normal, r.origin)) / denom;
  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  var intersection = ray_at(r, t);
  var planar_hitpt_vector = intersection - Q.xyz;
  var alpha = dot(w, cross(planar_hitpt_vector, v.xyz));
  var beta = dot(w, cross(u.xyz, planar_hitpt_vector));

  if (alpha < 0.0 || alpha > 1.0 || beta < 0.0 || beta > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (dot(normal, r.direction) > 0.0)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = intersection;
  record.normal = normal;
  record.hit_anything = true;
}

fn hit_triangle(r: ray, v0: vec3f, v1: vec3f, v2: vec3f, record: ptr<function, hit_record>, max: f32)
{
  var v1v0 = v1 - v0;
  var v2v0 = v2 - v0;
  var rov0 = r.origin - v0;

  var n = cross(v1v0, v2v0);
  var q = cross(rov0, r.direction);

  var d = 1.0 / dot(r.direction, n);

  var u = d * dot(-q, v2v0);
  var v = d * dot(q, v1v0);
  var t = d * dot(-n, rov0);

  if (u < 0.0 || u > 1.0 || v < 0.0 || (u + v) > 1.0)
  {
    record.hit_anything = false;
    return;
  }

  if (t < RAY_TMIN || t > max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = normalize(n);
  record.hit_anything = true;
}

fn hit_box(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, t_max: f32)
{
  var m = 1.0 / r.direction;
  var n = m * (r.origin - center);
  var k = abs(m) * rad;

  var t1 = -n - k;
  var t2 = -n + k;

  var tN = max(max(t1.x, t1.y), t1.z);
  var tF = min(min(t2.x, t2.y), t2.z);

  if (tN > tF || tF < 0.0)
  {
    record.hit_anything = false;
    return;
  }

  var t = tN;
  if (t < RAY_TMIN || t > t_max)
  {
    record.hit_anything = false;
    return;
  }

  record.t = t;
  record.p = ray_at(r, t);
  record.normal = -sign(r.direction) * step(t1.yzx, t1.xyz) * step(t1.zxy, t1.xyz);
  record.hit_anything = true;

  return;
}

fn hit_column(r: ray, center: vec3f, rad: vec3f, record: ptr<function, hit_record>, max: f32) {
  var r_d = vec3f(r.direction[0], 0.0, r.direction[2]);
  var ro_d = vec3f(r.origin[0], 0.0, r.origin[2]);
  var c_d = vec3f(center[0], 0.0, center[2]);
  var radius = f32(rad[1]);

  var a = dot(r_d, r_d);
  var b = 2*dot(r_d, (ro_d - c_d));
  var c = dot((ro_d - c_d), (ro_d - c_d)) - radius * radius;;
  var delta = pow(b,2) - 4*a*c;

  var t = f32(0.0);
  if (delta < 0) {
    record.hit_anything = false;
    return;
  }
  if (delta > 0) {
    var t1 = (-b + sqrt(delta))/(2*a);
    var t2 = (-b - sqrt(delta))/(2*a);
    if (t1 < t2 && t1 > RAY_TMIN) {
      t = t1;
    }
    else {
      t = t2;
    }
  }
  if (delta == 0) {
    t = -b/(2*a);
  }
  if (t < RAY_TMIN || t > max) {
      record.hit_anything = false;
      return;
  }

  var intersection = ray_at(r,t);
  var i_d = vec3f(intersection[0], 0.0, intersection[2]);
  var n = normalize(i_d - c_d);
  var normal = n;

  record.t = t;
  record.p = intersection;
  record.normal = normal;
  record.hit_anything = true;
}