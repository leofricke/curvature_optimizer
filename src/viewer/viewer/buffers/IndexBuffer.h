#pragma once

#include <pmp/viewers/gl.h>

class IndexBuffer
{
public:
    IndexBuffer();
    virtual ~IndexBuffer();

    IndexBuffer(const IndexBuffer&) = delete;
    IndexBuffer& operator=(const IndexBuffer&) = delete;

    IndexBuffer(IndexBuffer&& buffer) noexcept;
    IndexBuffer& operator=(IndexBuffer&& buffer) noexcept;

    void upload(const unsigned int* data, size_t size);

    void bind() const;
    void unbind() const;

private:
    GLuint m_handle;
    size_t m_size;
};