#ifndef TOOL_H
#define TOOL_H
#include<pch.h>
class Tool{
public:
    Tool(){}
    static inline float dotProduct(Vector3 a, Vector3 b){
        return a.dot(b);
    }
    static inline Vector3 crossProduct(Vector3 a, Vector3 b){
        return a.cross(b);
    }
    static inline float dotProduct(Vector2 a, Vector2 b){
        return a.dot(b);
    }
    static inline float normal(Vector2 a, Vector2 b){
        return a.dot(b);
    }


};
#endif // TOOL_H
