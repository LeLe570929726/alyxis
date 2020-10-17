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
    pub static mut IMPL_INS: &dyn MatImpl = &general::FPU_IMPL_INS;

    trait MatImpl {
        unsafe fn add_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]);
        unsafe fn add_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]);
        unsafe fn add_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]);
        unsafe fn sub_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]);
        unsafe fn sub_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]);
        unsafe fn sub_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]);
        unsafe fn mul_f32_mat2(&self, a: &mut [f32; 4], b: f32);
        unsafe fn mul_f32_mat3(&self, a: &mut [f32; 9], b: f32);
        unsafe fn mul_f32_mat4(&self, a: &mut [f32; 16], b: f32);
        unsafe fn mul_vec2_mat2(&self, a: &[f32; 4], b: &mut [f32; 4]);
        unsafe fn mul_vec3_mat3(&self, a: &[f32; 4], b: &mut [f32; 4]);
        unsafe fn mul_vec4_mat4(&self, a: &[f32; 4], b: &mut [f32; 4]);
        unsafe fn mul_mat2_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]);
        unsafe fn mul_mat3_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]);
        unsafe fn mul_mat4_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]);
        unsafe fn transpose_mat2(&self, a: &mut [f32; 4]);
        unsafe fn transpose_mat3(&self, a: &mut [f32; 9]);
        unsafe fn transpose_mat4(&self, a: &mut [f32; 16]);
        unsafe fn det_mat2(&self, a: &[f32; 4]) -> f32;
        unsafe fn det_mat3(&self, a: &[f32; 9]) -> f32;
        unsafe fn det_mat4(&self, a: &[f32; 16]) -> f32;
    }

    pub mod general {
        pub const FPU_IMPL_INS: FPUImpl = FPUImpl();

        pub struct FPUImpl();

        impl super::MatImpl for FPUImpl {
            unsafe fn add_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
                    *x += y;
                });
            }

            unsafe fn add_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]) {
                a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
                    *x += y;
                });
            }

            unsafe fn add_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]) {
                a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
                    *x += y;
                });
            }

            unsafe fn sub_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
                    *x -= y;
                });
            }

            unsafe fn sub_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]) {
                a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
                    *x -= y;
                });
            }

            unsafe fn sub_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]) {
                a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
                    *x -= y;
                });
            }

            unsafe fn mul_f32_mat2(&self, a: &mut [f32; 4], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x *= b;
                });
            }

            unsafe fn mul_f32_mat3(&self, a: &mut [f32; 9], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x *= b;
                });
            }

            unsafe fn mul_f32_mat4(&self, a: &mut [f32; 16], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x *= b;
                });
            }

            unsafe fn mul_vec2_mat2(&self, a: &[f32; 4], b: &mut [f32; 4]) {
                let src_vec = b.clone();
                for row in 0..2 {
                    let start_idx = row * 2;
                    let end_idx = start_idx + 2;
                    b[row] = a[start_idx..start_idx].iter().zip(src_vec.iter()).map(|(x, y)| {
                        x * y
                    }).sum();
                }
            }

            unsafe fn mul_vec3_mat3(&self, a: &[f32; 4], b: &mut [f32; 4]) {
                let src_vec = b.clone();
                for row in 0..3 {
                    let start_idx = row * 3;
                    let end_idx = start_idx + 3;
                    b[row] = a[start_idx..end_idx].iter().zip(src_vec.iter()).map(|(x, y)| {
                        x * y
                    }).sum();
                }
            }

            unsafe fn mul_vec4_mat4(&self, a: &[f32; 4], b: &mut [f32; 4]) {
                let src_vec = b.clone();
                for row in 0..4 {
                    let start_idx = row * 4;
                    let end_idx = start_idx + 4;
                    b[row] = a[start_idx..end_idx].iter().zip(src_vec.iter()).map(|(x, y)| {
                        x * y
                    }).sum();
                }
            }

            unsafe fn mul_mat2_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                let src_mat = a.clone();
                for row in 0..2 {
                    let start_idx = row * 2;
                    for col in 0..2 {
                        let idx = start_idx + col;
                        a[idx] = src_mat[start_idx] * b[col] +
                            src_mat[start_idx + 1] * b[col + 2];
                    }
                }
            }

            unsafe fn mul_mat3_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]) {
                let src_mat = a.clone();
                for row in 0..3 {
                    let start_idx = row * 3;
                    for col in 0..3 {
                        let idx = start_idx + col;
                        a[idx] = src_mat[start_idx] * b[col] +
                            src_mat[start_idx + 1] * b[col + 3] +
                            src_mat[start_idx + 2] * b[col + 6];
                    }
                }
            }

            unsafe fn mul_mat4_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]) {
                let src_mat = a.clone();
                for row in 0..4 {
                    let start_idx = row * 4;
                    for col in 0..4 {
                        let idx = start_idx + col;
                        a[idx] = src_mat[start_idx] * b[col] +
                            src_mat[start_idx + 1] * b[col + 4] +
                            src_mat[start_idx + 2] * b[col + 8] +
                            src_mat[start_idx + 3] * b[col + 12];
                    }
                }
            }

            unsafe fn transpose_mat2(&self, a: &mut [f32; 4]) {
                a[1] = a[2];
            }

            unsafe fn transpose_mat3(&self, a: &mut [f32; 9]) {
                a[1] = a[3];
                a[2] = a[6];
                a[5] = a[7];
            }

            unsafe fn transpose_mat4(&self, a: &mut [f32; 16]) {
                a[1] = a[4];
                a[2] = a[8];
                a[3] = a[12];
                a[6] = a[9];
                a[7] = a[13];
                a[11] = a[14];
            }

            unsafe fn det_mat2(&self, a: &[f32; 4]) -> f32 {
                a[0] * a[3] - a[1] * a[2]
            }

            unsafe fn det_mat3(&self, a: &[f32; 9]) -> f32 {
                let r1 = a[4] * a[8] - a[5] * a[7];
                let r2 = a[3] * a[8] - a[5] * a[6];
                let r3 = a[3] * a[7] - a[4] * a[6];
                a[0] * r1 - a[1] * r2 + a[2] * r3
            }

            unsafe fn det_mat4(&self, a: &[f32; 16]) -> f32 {
                let r11 = a[8] * a[13] - a[9] * a[12];
                let r12 = a[8] * a[14] - a[10] * a[12];
                let r13 = a[8] * a[15] - a[11] * a[12];
                let r14 = a[9] * a[14] - a[10] * a[13];
                let r15 = a[9] * a[15] - a[11] * a[13];
                let r16 = a[10] * a[15] - a[11] * a[14];
                let r21 = a[4] * r14 - a[5] * r12 + a[6] * r11;
                let r22 = a[4] * r15 - a[5] * r13 + a[7] * r11;
                let r23 = a[4] * r16 - a[6] * r13 + a[7] * r12;
                let r24 = a[5] * r16 - a[6] * r15 + a[7] * r14;
                a[0] * r24 - a[1] * r23 + a[2] * r22 - a[3] * r21
            }
        }
    }
}
