#pragma once

#include "VertexBuffer.h"
#include <pmp/viewers/gl.h>
#include <vector>

class VertexArray
{
public:
    explicit VertexArray();
    virtual ~VertexArray();

    VertexArray(const VertexArray&) = delete;
    VertexArray& operator=(const VertexArray&) = delete;

    VertexArray(VertexArray&& array) noexcept;
    VertexArray& operator=(VertexArray&& array) noexcept;

    void addBuffer(GLint size, GLuint divisor, GLenum glType);

    const VertexBuffer& getBuffer(size_t index) const;
    VertexBuffer& getBuffer(size_t index);
    size_t numBuffers() const;

    void bind() const;
    void unbind() const;

private:
    GLuint m_handle;
    std::vector<VertexBuffer> m_buffers;
};