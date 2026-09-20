#include "VertexArray.h"

VertexArray::VertexArray()
    : m_handle(0)
{
}

VertexArray::~VertexArray()
{
    if (m_handle)
    {
        glDeleteVertexArrays(1, &m_handle);
    }
}

VertexArray::VertexArray(VertexArray&& array) noexcept
{
    m_handle = array.m_handle;
    m_buffers = std::move(array.m_buffers);
    array.m_handle = 0;
}

VertexArray& VertexArray::operator=(VertexArray&& array) noexcept
{
    m_handle = array.m_handle;
    m_buffers = std::move(array.m_buffers);
    array.m_handle = 0;
    return *this;
}

void VertexArray::addBuffer(GLint size, GLuint divisor, GLenum glType)
{
    if (!m_handle)
    {
        glGenVertexArrays(1, &m_handle);
    }

    const VertexAttribInfo attribInfo = { .m_index = static_cast<GLuint>(m_buffers.size()), .m_size = size, .m_divisor = divisor, .m_glType = glType };
    m_buffers.emplace_back(attribInfo);
}

const VertexBuffer& VertexArray::getBuffer(size_t index) const
{
    return m_buffers[index];
}

VertexBuffer& VertexArray::getBuffer(size_t index)
{
    return m_buffers[index];
}

size_t VertexArray::numBuffers() const
{
    return m_buffers.size();
}

void VertexArray::bind() const
{
    if (m_handle)
    {
        glBindVertexArray(m_handle);
    }
}

void VertexArray::unbind() const
{
    if (m_handle)
    {
        glBindVertexArray(0);
    }
}
