#include "VertexBuffer.h"

VertexBuffer::VertexBuffer(const VertexAttribInfo& attribInfo)
    : m_handle(0), m_size(0), m_attribInfo(attribInfo)
{
}

VertexBuffer::~VertexBuffer()
{
    if (m_handle)
    {
        glDeleteBuffers(1, &m_handle);
    }
}

VertexBuffer::VertexBuffer(VertexBuffer&& buffer) noexcept
{
    m_handle = buffer.m_handle;
    m_size = buffer.m_size;
    m_attribInfo = buffer.m_attribInfo;
    buffer.m_handle = 0;
}

VertexBuffer& VertexBuffer::operator=(VertexBuffer&& buffer) noexcept
{
    m_handle = buffer.m_handle;
    m_size = buffer.m_size;
    m_attribInfo = buffer.m_attribInfo;
    buffer.m_handle = 0;
    return *this;
}

void VertexBuffer::upload(const void* data, const size_t size)
{
    if (!m_handle)
    {
        glGenBuffers(1, &m_handle);
    }

    bind();
    m_size = size;
    glBufferData(GL_ARRAY_BUFFER, size, data, GL_STATIC_DRAW);
    glEnableVertexAttribArray(m_attribInfo.m_index);
    glVertexAttribPointer(m_attribInfo.m_index, m_attribInfo.m_size, m_attribInfo.m_glType, GL_FALSE, 0, nullptr);
    if (m_attribInfo.m_divisor > 0)
    {
        glVertexAttribDivisor(m_attribInfo.m_index, m_attribInfo.m_divisor);
    }
    unbind();
}

void VertexBuffer::bind() const
{
    if (m_handle)
    {
        glBindBuffer(GL_ARRAY_BUFFER, m_handle);
    }
}

void VertexBuffer::unbind() const
{
    if (m_handle)
    {
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
}