#include <poplar/Vertex.hpp>
#include <poplar/FieldTypes.hpp>
#include <poplar/HalfFloat.hpp>
#include <poplar/StackSizeDefs.hpp>
#include <poplar/AvailableVTypes.h>
#include "poplar/TileConstants.hpp"
#include <print.h>
#include <type_traits>


class AddMod : public poplar::Vertex {
public:
    int mod;
    poplar::InOut<int> val;
    bool compute() {
            *val = (*val + 1) % mod;
            return true;
    }
};