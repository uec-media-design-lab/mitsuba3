#pragma once

#include <math.h>
#include <mitsuba/render/optix/common.h>
#include <mitsuba/render/optix/math.cuh>

struct OptixMyHalfSphereData {
    optix::BoundingBox3f bbox;
    optix::Transform4f to_object;
    optix::Vector3f center;
    float radius;
};

#ifdef __CUDACC__

extern "C" __global__ void __intersection__myhalfsphere() {
    const OptixHitGroupData *sbt_data = (OptixHitGroupData*) optixGetSbtDataPointer();
    OptixMyHalfSphereData *sphere = (OptixMyHalfSphereData*) sbt_data->data;

    // Ray in instance-space
    Ray3f ray = get_ray();

    /* 以下追加・変更 */
    Ray3f ray_local = sphere->to_object.transform_ray(ray);

    // We define a plane which is perpendicular to the ray direction and
    // contains the sphere center and intersect it. We then solve the ray-sphere
    // intersection as if the ray origin was this new intersection point. This
    // additional step makes the whole intersection routine numerically more
    // robust.

    Vector3f l = ray.o - sphere->center;
    Vector3f d = ray.d;
    float plane_t = dot(-l, d) / norm(d);

    // Ray is perpendicular to plane
    if (plane_t == 0 && ray.o != sphere->center)
        return;

    Vector3f plane_p = ray(plane_t);

    // Intersection with plane outside of the sphere
    if (norm(plane_p - sphere->center) > sphere->radius)
        return;

    Vector3f o = plane_p - sphere->center;

    float A = squared_norm(d);
    float B = 2.f * dot(o, d);
    float C = squared_norm(o) - sqr(sphere->radius);

    float near_t, far_t;
    bool solution_found = solve_quadratic(A, B, C, near_t, far_t);

    near_t += plane_t;
    far_t += plane_t;

    // Sphere doesn't intersect with the segment on the ray
    bool out_bounds = !(near_t <= ray.maxt && far_t >= 0.f); // NaN-aware conditionals

    // Sphere fully contains the segment of the ray
    bool in_bounds = near_t < 0.f && far_t > ray.maxt;

    /* 以下追加・変更 */
    float z_pos_near = ray_local.o.z() + ray_local.d.z() * near_t,
          z_pos_far  = ray_local.o.z() + ray_local.d.z() * far_t;

    float t = (near_t < 0.f || z_pos_near < 0.f ? far_t: near_t);

    if (solution_found && !out_bounds && !in_bounds && (z_pos_near >= 0.f || z_pos_far >= 0.f))
        optixReportIntersection(t, OPTIX_HIT_KIND_TRIANGLE_FRONT_FACE);
}

extern "C" __global__ void __closesthit__myhalfsphere() {
    const OptixHitGroupData *sbt_data = (OptixHitGroupData *) optixGetSbtDataPointer();
    set_preliminary_intersection_to_payload(
        optixGetRayTmax(), Vector2f(), 0, sbt_data->shape_registry_id);
}
#endif
