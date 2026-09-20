#pragma once

#include <pmp/viewers/gl.h>

struct VertexAttribInfo
{
    GLuint m_index;
    GLint m_size;
    GLuint m_divisor;
    GLenum m_glType;

    size_t sizeOfGLType() const
    {
        switch (m_glType)
        {
        case GL_FLOAT:
            return 4;
        case GL_UNSIGNED_INT:
            return 4;
        case GL_UNSIGNED_BYTE:
            return 1;
        }
        return 0;
    }
};

class VertexBuffer
{
public:
    explicit VertexBuffer(const VertexAttribInfo& attribInfo);
    virtual ~VertexBuffer();

    VertexBuffer(const VertexBuffer&) = delete;
    VertexBuffer& operator=(const VertexBuffer&) = delete;

    VertexBuffer(VertexBuffer&& buffer) noexcept;
    VertexBuffer& operator=(VertexBuffer&& buffer) noexcept;

    void upload(const void* data, size_t size);

    void bind() const;
    void unbind() const;

private:
    GLuint m_handle;
    size_t m_size;
    VertexAttribInfo m_attribInfo;
};