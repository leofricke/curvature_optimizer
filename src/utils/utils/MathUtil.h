#pragma once

#include <pmp/surface_mesh.h>

pmp::dmat2 rotationMatrix(double angle);
pmp::dmat4 rotationMatrix(const pmp::dvec3& axis, double angle);
pmp::dvec3 affineTransform(const pmp::dmat4& m, const pmp::dvec3& v);
pmp::dvec3 linearTransform(const pmp::dmat4& m, const pmp::dvec3& v);

pmp::dmat2 inverse(const pmp::dmat2& m);
