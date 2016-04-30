#include "scene/3d/skeleton.h"
#include "scene/3d/camera.h"

#define GME_MATH_TAU 3.14159265358979323846 * 2

class GodotMathExtension : public Object {

	OBJ_TYPE(GodotMathExtension, Object);

	static GodotMathExtension *singleton;
public:
	static _FORCE_INLINE_ Vector3 spherical_to_local_position(real_t p_theta, real_t p_phi) {
		Vector3 res;

		real_t sin_theta = sin(p_theta);
		real_t cos_theta = cos(p_theta);
		real_t sin_phi = sin(p_phi);
		real_t cos_phi = cos(p_phi);

		res.z = sin_theta * cos_phi;
		res.y = cos_theta;
		res.x = sin_theta * sin_phi;

		return res;
	};

	static _FORCE_INLINE_ Quat quat_from_radians(Vector3 p_radians) {
		real_t pitch_radians = p_radians.x * 0.5;
		real_t yaw_radians = p_radians.y * 0.5;
		real_t roll_radians = p_radians.z * 0.5;

		real_t sin_pitch = sin(pitch_radians);
		real_t cos_pitch = cos(pitch_radians);

		real_t sin_yaw = sin(yaw_radians);
		real_t cos_yaw = cos(yaw_radians);

		real_t sin_roll = sin(roll_radians);
		real_t cos_roll = cos(roll_radians);

		return Quat(sin_yaw * cos_pitch * sin_roll + cos_yaw * sin_pitch * cos_roll, sin_yaw * cos_pitch * cos_roll - cos_yaw * sin_pitch * sin_roll, cos_yaw * cos_pitch * sin_roll - sin_yaw * sin_pitch * cos_roll, cos_yaw * cos_pitch * cos_roll + sin_yaw * sin_pitch * sin_roll);
	};

	static _FORCE_INLINE_ real_t ease_in(real_t t) {
		return sin(t * Math_PI * 0.5);
	}

	static _FORCE_INLINE_ real_t ease_out(real_t t) {
		return cos(t * Math_PI * 0.5);
	}

	static _FORCE_INLINE_ real_t exponetial(real_t t) {
		return t * t;
	}

	static _FORCE_INLINE_ real_t smooth_step(real_t t) {
		return t * t * (3.0 - 2.0 * t);
	}

	static _FORCE_INLINE_ real_t smoother_step(real_t t) {
		return t * t * t * (t * (6.0 * t - 15.0) + 10.0);
	}

	static _FORCE_INLINE_ real_t camera_get_position_distance(const Camera *p_camera, const Vector3 &p_pos) {
		Transform t = p_camera->get_global_transform();
		Vector3 axis = -Vector3(t.basis[2].x, t.basis[2].y, t.basis[2].z);
		Vector3 eyedir = axis.normalized();
		return eyedir.dot(p_pos) - (eyedir.dot(t.origin));
	}

	static _FORCE_INLINE_ Vector2 get_2d_position_from_3d_position_with_screen_limits(const Camera *p_camera, const Vector3 &p_position_3d,
		const Vector2 &screen_size, const Vector2 &screen_center,
		const Vector2 &screen_mins, const Vector2 &screen_max) {

		bool is_behind = camera_get_position_distance(p_camera, p_position_3d) < 0;
		
		Vector2 screen_pos = p_camera->unproject_position(p_position_3d);

		Vector2 screen_bounds_min = screen_center - screen_mins;
		Vector2 screen_bounds_max = screen_center - (screen_size - screen_max);

		if (is_behind == false &&
			(screen_pos.x>(screen_mins.x) && screen_pos.x < (screen_max.x)) &&
			(screen_pos.y>(screen_mins.y) && screen_pos.y < (screen_max.y))) {
			return Vector2();
		}
		else {
			int rotation = 270;

			if (is_behind == true)
				rotation = 90;
			else
				rotation = 270;

			screen_pos = screen_pos - screen_center;

			real_t angle = atan2(screen_pos.y, screen_pos.x);
			angle = angle - rotation * ((Math_PI * 2) / 360);

			real_t angle_cos = cos(angle);
			real_t angle_sin = sin(angle);

			real_t m = angle_cos / angle_sin;
			
			if (angle_cos > 0)
				screen_pos = Vector2(screen_bounds_min.y / m, screen_bounds_min.y);
			else
				screen_pos = Vector2(-screen_bounds_max.y / m, -screen_bounds_max.y);
		
			if (screen_pos.x > screen_bounds_max.x)
				screen_pos = Vector2(screen_bounds_max.x, screen_bounds_max.x*m);
			else if(screen_pos.x < -screen_bounds_min.x)
				screen_pos = Vector2(-screen_bounds_min.x, -screen_bounds_min.x*m);

			screen_pos = screen_pos + screen_center;
			screen_pos.y = screen_size.y - screen_pos.y;
		}
		return screen_pos;
	};

	static _FORCE_INLINE_ Vector2 get_2d_position_from_3d_position(const Camera *p_camera, const Vector3 &p_position_3d) {
		Vector2 screen_pos = p_camera->unproject_position(p_position_3d);
		return screen_pos;
	}

	static _FORCE_INLINE_ real_t clamp_angle(real_t val, real_t ang_min, real_t ang_max) {
		if (val < -360)
			val += 360;
		if (val > 360)
			val -= 360;
		return CLAMP(val, ang_min, ang_max);
	}

	static _FORCE_INLINE_ Vector3 adjust_facing(const Vector3 &p_facing, const Vector3 &p_target, const real_t &p_step, const real_t &p_adjust_rate, const Vector3& p_current_gn) {
		Vector3 n = p_target; //normal
		Vector3 t = n.cross(p_current_gn).normalized();

		real_t x = n.dot(p_facing);
		real_t y = t.dot(p_facing);

		real_t ang = atan2(y, x);

		if (abs(ang) < 0.001)
			return p_facing;

		real_t s = ang < 0.0 ? -1.0 : (ang > 0.0 ? +1.0 : 0.0); // Sign
		ang = ang * s;
		real_t turn = ang * p_adjust_rate * p_step;
		real_t a;
		
		if (ang < turn) {
			a = ang;
		}
		else {
			a = turn;
		}

		ang = (ang - a) * s;

		return ((n * cos(ang)) + (t * sin(ang))) * p_facing.length();
	}

	static _FORCE_INLINE_ Transform rotate_around(Transform p_transform, Vector3 p_point, Vector3 p_axis, real_t p_angle) {
		Vector3 vector = p_transform.origin;
		Quat rotation = Quat(p_axis, p_angle);
		Vector3 vector2 = vector - p_point;
		Quat quat = rotation * vector2;
		vector = p_point + vector2;
		p_transform.origin = vector;

		return p_transform.rotated(p_axis, p_angle * 0.0174532924f);
	}
};

class _GodotMathExtension : public Object {

	OBJ_TYPE(_GodotMathExtension, Object);

	static _GodotMathExtension *singleton;
protected:
	static void _bind_methods();

public:
	static _GodotMathExtension *get_singleton();

	Vector3 spherical_to_local_position(real_t p_theta, real_t p_phi);
	Quat quat_from_radians(Vector3 p_radians);
	real_t ease_in(real_t t);
	real_t ease_out(real_t t);
	real_t exponetial(real_t t);
	real_t smooth_step(real_t t);
	real_t smoother_step(real_t t);
	real_t camera_get_position_distance(const Object *p_camera, const Vector3 &p_pos);

	Vector2 get_2d_position_from_3d_position_with_screen_limits(const Object *p_camera, const Vector3 &p_position_3d,
	const Vector2 &screen_size, const Vector2 &screen_center,
	const Vector2 &screen_mins, const Vector2 &screen_max);

	Vector2 get_2d_position_from_3d_position(const Object *p_camera, const Vector3 &p_position_3d);
	real_t clamp_angle(real_t val, real_t ang_min, real_t ang_max);
	Vector3 adjust_facing(const Vector3 &p_facing, const Vector3 &p_target, const real_t &p_step, const real_t &p_adjust_rate, const Vector3& p_current_gn);
	Transform rotate_around(Transform p_transform, Vector3 p_point, Vector3 p_axis, real_t p_angle);

	_GodotMathExtension();
};

static _GodotMathExtension *_godot_math_extension = NULL;