#pragma once

#include <GL/glew.h>
#include <GL/GLU.h>

GLuint CompileGLShader(const char* pchShaderName, const char* pchVertexShader, const char* pchFragmentShader);

void dprintf(const char* fmt, ...);

static bool g_bPrintf = true;
