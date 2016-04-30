#include "register_types.h"

#include "globals.h"
#include "godot_math_extension.h"

void register_godot_math_extension_types() {
	_godot_math_extension = memnew(_GodotMathExtension);

	Globals::get_singleton()->add_singleton(Globals::Singleton("GodotMathExtension", _GodotMathExtension::get_singleton()));
}
void unregister_godot_math_extension_types() {
	memdelete(_godot_math_extension);
}
