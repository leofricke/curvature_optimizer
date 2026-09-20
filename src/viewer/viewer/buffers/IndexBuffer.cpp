#include "IndexBuffer.h"

IndexBuffer::IndexBuffer()
    : m_handle(0), m_size(0)
{
}

IndexBuffer::~IndexBuffer()
{
    if (m_handle)
    {
        glDeleteBuffers(1, &m_handle);
    }
}
IndexBuffer::IndexBuffer(IndexBuffer&& buffer) noexcept
{
    m_handle = buffer.m_handle;
    m_size = buffer.m_size;
    buffer.m_handle = 0;
}

IndexBuffer& IndexBuffer::operator=(IndexBuffer&& buffer) noexcept
{
    m_handle = buffer.m_handle;
    m_size = buffer.m_size;
    buffer.m_handle = 0;
    return *this;
}

void IndexBuffer::upload(const unsigned int* data, const size_t size)
{
    if (!m_handle)
    {
        glGenBuffers(1, &m_handle);
    }

    bind();
    m_size = size;
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, size * sizeof(unsigned int), data, GL_STATIC_DRAW);
    unbind();
}

void IndexBuffer::bind() const
{
    if (m_handle)
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_handle);
    }
}

void IndexBuffer::unbind() const
{
    if (m_handle)
    {
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
}