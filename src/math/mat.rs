#[repr(align(32))]
#[derive(Clone)]
pub struct Mat2([f32; 4]);

#[repr(align(32))]
#[derive(Clone)]
pub struct Mat3([f32; 9]);

#[repr(align(32))]
#[derive(Clone)]
pub struct Mat4([f32; 16]);

impl Mat2 {
    pub fn new(m11: f32, m12: f32, m21: f32, m22: f32) -> Self {
        Self([m11, m12, m21, m22])
    }
}

impl Mat3 {
    pub fn new(m11: f32, m12: f32, m13: f32, m21: f32, m22: f32, m23: f32, m31: f32, m32: f32, m33: f32) -> Self {
        Self([m11, m12, m13, m21, m22, m23, m31, m32, m33])
    }
}

impl Mat4 {
    pub fn new(m11: f32, m12: f32, m13: f32, m14: f32, m21: f32, m22: f32, m23: f32, m24: f32, m31: f32, m32: f32, m33: f32, m34: f32, m41: f32, m42: f32, m43: f32, m44: f32) -> Self {
        Self([m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44])
    }
}

mod impls {
    trait MatImpl {

    }
}