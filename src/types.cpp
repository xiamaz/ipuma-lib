#include "types.h"
#include <string>

namespace ipu {
std::string vertexTypeToString(VertexType v) { return typeString[static_cast<int>(v)]; }
}