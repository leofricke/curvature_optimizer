static const char* point_vshader =
#ifndef __EMSCRIPTEN__
    "#version 330"
#else
    "#version 300 es"
#endif
    R"glsl(
layout (location=0) in vec3 v_position;
layout (location=1) in vec3 v_point;
layout (location=2) in vec3 v_color;

out vec3 u_color;

uniform float u_screenWidth;
uniform float u_screenHeight;
uniform float u_pointSize;
uniform mat4 u_mvpMatrix;

void main()
{
    vec2 resolution = vec2(u_screenWidth, u_screenHeight);

    vec4 pos = u_mvpMatrix * vec4(v_point, 1.0);
    vec2 screen = resolution * (0.5 * pos.xy / pos.w + 0.5);

    vec2 p = screen + u_pointSize * v_position.xy;

    u_color = v_color;
    gl_Position = vec4(pos.w * ((2.0 * p) / resolution - 1.0), pos.z, pos.w);
}
)glsl";

static const char* point_fshader =
#ifndef __EMSCRIPTEN__
    "#version 330\n"
#else
    "#version 300 es\n"
#endif
    R"glsl(
precision mediump float;

in vec3 u_color;

out vec4 f_color;

void main()
{
    f_color = vec4(u_color, 1.0);
}
)glsl";
