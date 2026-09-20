static const char* line_vshader =
#ifndef __EMSCRIPTEN__
    "#version 330"
#else
    "#version 300 es"
#endif
    R"glsl(
layout (location=0) in vec3 v_position;
layout (location=1) in vec3 v_start;
layout (location=2) in vec3 v_end;
layout (location=3) in vec3 v_color;

out vec3 u_color;

uniform float u_screenWidth;
uniform float u_screenHeight;
uniform float u_lineWidth;
uniform mat4 u_mvpMatrix;

void main()
{
    vec2 resolution = vec2(u_screenWidth, u_screenHeight);

    vec4 clipStart = u_mvpMatrix * vec4(v_start, 1.0);
    vec4 clipEnd = u_mvpMatrix * vec4(v_end, 1.0);
    
    vec2 screenStart = resolution * (0.5 * clipStart.xy / clipStart.w + 0.5);
    vec2 screenEnd = resolution * (0.5 * clipEnd.xy / clipEnd.w + 0.5);

    vec2 dir = normalize(screenEnd - screenStart);
    vec2 dir2 = vec2(-dir.y, dir.x);

    vec2 pStart = screenStart + u_lineWidth * (v_position.x * dir + v_position.y * dir2);
    vec2 pEnd = screenEnd + u_lineWidth * (v_position.x * dir + v_position.y * dir2);

    vec4 clip = mix(clipStart, clipEnd, v_position.z);
    vec2 p = mix(pStart, pEnd, v_position.z);

    u_color = v_color;
    gl_Position = vec4(clip.w * ((2.0 * p) / resolution - 1.0), clip.z, clip.w);
}
)glsl";

static const char* line_fshader =
#ifndef __EMSCRIPTEN__
    "#version 330\n"
#else
    "#version 300 es\n"
#endif
    R"glsl(
precision mediump float;

in vec3 u_color;

out vec4 f_color;

uniform float u_lineAlpha;

void main()
{
    f_color = vec4(u_color, u_lineAlpha);
}
)glsl";
