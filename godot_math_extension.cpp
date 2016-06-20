#include "godot_math_extension.h"

_GodotMathExtension *_GodotMathExtension::singleton = NULL;

void _GodotMathExtension::_bind_methods() {
	ObjectTypeDB::bind_method(_MD("spherical_to_local_position:Vector3", "theta", "phi"), &_GodotMathExtension::spherical_to_local_position);
	ObjectTypeDB::bind_method(_MD("quat_from_radians:Quat", "radians"), &_GodotMathExtension::quat_from_radians);

	ObjectTypeDB::bind_method(_MD("ease_in", "t"), &_GodotMathExtension::ease_in);
	ObjectTypeDB::bind_method(_MD("ease_out", "t"), &_GodotMathExtension::ease_out);
	ObjectTypeDB::bind_method(_MD("exponetial", "t"), &_GodotMathExtension::exponetial);
	ObjectTypeDB::bind_method(_MD("smooth_step", "t"), &_GodotMathExtension::smooth_step);
	ObjectTypeDB::bind_method(_MD("smoother_step", "t"), &_GodotMathExtension::smoother_step);

	ObjectTypeDB::bind_method(_MD("camera_get_position_distance", "camera", "pos"), &_GodotMathExtension::camera_get_position_distance);
	//ObjectTypeDB::bind_method(_MD("get_2d_position_from_3d_position_with_screen_limits:Vector2", "camera", "position_3d", "screen_size", "screen_center", "screen_mins", "screen_max"), &_GodotMathExtension::get_2d_position_from_3d_position_with_screen_limits);
	ObjectTypeDB::bind_method(_MD("get_2d_position_from_3d_position:Vector2", "camera", "position_3d"), &_GodotMathExtension::get_2d_position_from_3d_position);
	ObjectTypeDB::bind_method(_MD("clamp_angle", "val", "ang_min", "ang_max"), &_GodotMathExtension::clamp_angle);
	ObjectTypeDB::bind_method(_MD("adjust_facing:Vector3", "facing", "target", "step", "adjust_rate", "current_gn"), &_GodotMathExtension::adjust_facing);
	ObjectTypeDB::bind_method(_MD("rotate_around:Transform", "transform", "point", "axis", "angle"), &_GodotMathExtension::rotate_around);
	ObjectTypeDB::bind_method(_MD("base_log:float", "float", "float"), &_GodotMathExtension::base_log);

	BIND_CONSTANT(GME_MATH_TAU);
}

_GodotMathExtension *_GodotMathExtension::get_singleton() {
	return singleton;
}

Vector3 _GodotMathExtension::spherical_to_local_position(real_t p_theta, real_t p_phi) {
	return GodotMathExtension::spherical_to_local_position(p_theta, p_phi);
}

Quat _GodotMathExtension::quat_from_radians(Vector3 p_radians) {
	return GodotMathExtension::quat_from_radians(p_radians);
}

real_t _GodotMathExtension::ease_in(real_t t) {
	return GodotMathExtension::ease_in(t);
}

real_t _GodotMathExtension::ease_out(real_t t) {
	return GodotMathExtension::ease_out(t);
}

real_t _GodotMathExtension::exponetial(real_t t) {
	return GodotMathExtension::exponetial(t);
}

real_t _GodotMathExtension::smooth_step(real_t t) {
	return GodotMathExtension::smooth_step(t);
}

real_t _GodotMathExtension::smoother_step(real_t t) {
	return GodotMathExtension::smoother_step(t);
}

real_t _GodotMathExtension::camera_get_position_distance(const Object *p_camera, const Vector3 &p_pos) {
	if (p_camera) {
		const Camera *camera = p_camera->cast_to<Camera>();
		if (camera)
			return GodotMathExtension::camera_get_position_distance(camera, p_pos);
	}

	return 0.0f;
}

Vector2 _GodotMathExtension::get_2d_position_from_3d_position_with_screen_limits(const Object *p_camera, const Vector3 &p_position_3d,
	const Vector2 &screen_size, const Vector2 &screen_center,
	const Vector2 &screen_mins, const Vector2 &screen_max) {

	if (p_camera) {
		const Camera *camera = p_camera->cast_to<Camera>();
		if (camera)
			return GodotMathExtension::get_2d_position_from_3d_position_with_screen_limits(camera, p_position_3d, screen_size, screen_center, screen_mins, screen_max);
	}

	return Vector2();
}

Vector2 _GodotMathExtension::get_2d_position_from_3d_position(const Object *p_camera, const Vector3 &p_position_3d) {
	if (p_camera) {
		const Camera *camera = p_camera->cast_to<Camera>();
		if (camera)
			return GodotMathExtension::get_2d_position_from_3d_position(camera, p_position_3d);
	}

	return Vector2();
}

real_t _GodotMathExtension::clamp_angle(real_t val, real_t ang_min, real_t ang_max) {
	return GodotMathExtension::clamp_angle(val, ang_min, ang_max);
}

Vector3 _GodotMathExtension::adjust_facing(const Vector3 &p_facing, const Vector3 &p_target, const real_t &p_step, const real_t &p_adjust_rate, const Vector3& p_current_gn) {
	return GodotMathExtension::adjust_facing(p_facing, p_target, p_step, p_adjust_rate, p_current_gn);
}

Transform _GodotMathExtension::rotate_around(Transform p_transform, Vector3 p_point, Vector3 p_axis, real_t p_angle) {
	return GodotMathExtension::rotate_around(p_transform, p_point, p_axis, p_angle);
};

float _GodotMathExtension::base_log(float a, float new_base) {
	if (new_base == 1.0) {
		return NAN;
	}
	if (a != 1.0 && (new_base == 0.0 || Math::is_inf(new_base))) {
		return NAN;
	}
		
	return Math::log(a) / Math::log(new_base);
}

_GodotMathExtension::_GodotMathExtension() {
	singleton = this;
}