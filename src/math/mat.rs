use super::*;

#[repr(align(32))]
#[derive(Clone)]
pub struct Mat2(pub(super) [f32; 4]);

#[repr(align(32))]
#[derive(Clone)]
pub struct Mat3(pub(super) [f32; 9]);

#[repr(align(32))]
#[derive(Clone)]
pub struct Mat4(pub(super) [f32; 16]);

pub trait SquareMat {
    type Vec;

    fn add(self, other: &Self) -> Self;
    fn sub(self, other: &Self) -> Self;
    fn mul_f32(self, other: f32) -> Self;
    fn mul_vec(&self, other: Self::Vec) -> Self::Vec;
    fn mul_mat(self, other: &Self) -> Self;
    fn div(self, other: f32) -> Self;
    fn div_fast(self, other: f32) -> Self;
    fn transpose(self) -> Self;
    fn det(&self) -> f32;
}

impl Mat2 {
    pub fn new(m11: f32, m12: f32, m21: f32, m22: f32) -> Self {
        Self([m11, m12, m21, m22])
    }
}

impl SquareMat for Mat2 {
    type Vec = vec::Vec2;

    #[inline]
    fn add(mut self, other: &Self) -> Self {
        unsafe {
            impls::IMPL_INS.add_mat2(&mut self.0, &other.0);
        }
        self
    }

    #[inline]
    fn sub(mut self, other: &Self) -> Self {
        unsafe {
            impls::IMPL_INS.sub_mat2(&mut self.0, &other.0);
        }
        self
    }

    #[inline]
    fn mul_f32(mut self, other: f32) -> Self {
        unsafe {
            impls::IMPL_INS.mul_f32_mat2(&mut self.0, other);
        }
        self
    }

    #[inline]
    fn mul_vec(&self, mut other: Self::Vec) -> Self::Vec {
        unsafe {
            impls::IMPL_INS.mul_vec2_mat2(&self.0, &mut other.0);
        }
        other
    }

    #[inline]
    fn mul_mat(mut self, other: &Self) -> Self {
        unsafe {
            impls::IMPL_INS.mul_mat2_mat2(&mut self.0, &other.0);
        }
        self
    }

    #[inline]
    fn div(mut self, other: f32) -> Self {
        unsafe {
            impls::IMPL_INS.div_mat2(&mut self.0, other);
        }
        self
    }

    #[inline]
    fn div_fast(mut self, other: f32) -> Self {
        unsafe {
            impls::IMPL_INS.div_fast_mat2(&mut self.0, other);
        }
        self
    }

    #[inline]
    fn transpose(mut self) -> Self {
        unsafe {
            impls::IMPL_INS.transpose_mat2(&mut self.0);
        }
        self
    }

    #[inline]
    fn det(&self) -> f32 {
        unsafe {
            impls::IMPL_INS.det_mat2(&self.0)
        }
    }
}

impl Mat3 {
    pub fn new(m11: f32, m12: f32, m13: f32, m21: f32, m22: f32, m23: f32, m31: f32, m32: f32, m33: f32) -> Self {
        Self([m11, m12, m13, m21, m22, m23, m31, m32, m33])
    }
}

impl SquareMat for Mat3 {
    type Vec = vec::Vec3;

    #[inline]
    fn add(mut self, other: &Self) -> Self {
        unsafe {
            impls::IMPL_INS.add_mat3(&mut self.0, &other.0);
        }
        self
    }

    #[inline]
    fn sub(mut self, other: &Self) -> Self {
        unsafe {
            impls::IMPL_INS.sub_mat3(&mut self.0, &other.0);
        }
        self
    }

    #[inline]
    fn mul_f32(mut self, other: f32) -> Self {
        unsafe {
            impls::IMPL_INS.mul_f32_mat3(&mut self.0, other);
        }
        self
    }

    #[inline]
    fn mul_vec(&self, mut other: Self::Vec) -> Self::Vec {
        unsafe {
            impls::IMPL_INS.mul_vec3_mat3(&self.0, &mut other.0);
        }
        other
    }

    #[inline]
    fn mul_mat(mut self, other: &Self) -> Self {
        unsafe {
            impls::IMPL_INS.mul_mat3_mat3(&mut self.0, &other.0);
        }
        self
    }

    #[inline]
    fn div(mut self, other: f32) -> Self {
        unsafe {
            impls::IMPL_INS.div_mat3(&mut self.0, other);
        }
        self
    }

    #[inline]
    fn div_fast(mut self, other: f32) -> Self {
        unsafe {
            impls::IMPL_INS.div_fast_mat3(&mut self.0, other);
        }
        self
    }

    #[inline]
    fn transpose(mut self) -> Self {
        unsafe {
            impls::IMPL_INS.transpose_mat3(&mut self.0);
        }
        self
    }

    #[inline]
    fn det(&self) -> f32 {
        unsafe {
            impls::IMPL_INS.det_mat3(&self.0)
        }
    }
}

impl Mat4 {
    pub fn new(m11: f32, m12: f32, m13: f32, m14: f32, m21: f32, m22: f32, m23: f32, m24: f32, m31: f32, m32: f32, m33: f32, m34: f32, m41: f32, m42: f32, m43: f32, m44: f32) -> Self {
        Self([m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44])
    }
}

impl SquareMat for Mat4 {
    type Vec = vec::Vec4;

    #[inline]
    fn add(mut self, other: &Self) -> Self {
        unsafe {
            impls::IMPL_INS.add_mat4(&mut self.0, &other.0);
        }
        self
    }

    #[inline]
    fn sub(mut self, other: &Self) -> Self {
        unsafe {
            impls::IMPL_INS.sub_mat4(&mut self.0, &other.0);
        }
        self
    }

    #[inline]
    fn mul_f32(mut self, other: f32) -> Self {
        unsafe {
            impls::IMPL_INS.mul_f32_mat4(&mut self.0, other);
        }
        self
    }

    #[inline]
    fn mul_vec(&self, mut other: Self::Vec) -> Self::Vec {
        unsafe {
            impls::IMPL_INS.mul_vec4_mat4(&self.0, &mut other.0);
        }
        other
    }

    #[inline]
    fn mul_mat(mut self, other: &Self) -> Self {
        unsafe {
            impls::IMPL_INS.mul_mat4_mat4(&mut self.0, &other.0);
        }
        self
    }

    #[inline]
    fn div(mut self, other: f32) -> Self {
        unsafe {
            impls::IMPL_INS.div_mat4(&mut self.0, other);
        }
        self
    }

    #[inline]
    fn div_fast(mut self, other: f32) -> Self {
        unsafe {
            impls::IMPL_INS.div_fast_mat4(&mut self.0, other);
        }
        self
    }

    #[inline]
    fn transpose(mut self) -> Self {
        unsafe {
            impls::IMPL_INS.transpose_mat4(&mut self.0);
        }
        self
    }

    #[inline]
    fn det(&self) -> f32 {
        unsafe {
            impls::IMPL_INS.det_mat4(&self.0)
        }
    }
}

mod impls {
    pub static mut IMPL_INS: &dyn MatImpl = &general::FPU_IMPL_INS;

    pub trait MatImpl {
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
        unsafe fn mul_vec3_mat3(&self, a: &[f32; 9], b: &mut [f32; 4]);
        unsafe fn mul_vec4_mat4(&self, a: &[f32; 16], b: &mut [f32; 4]);
        unsafe fn mul_mat2_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]);
        unsafe fn mul_mat3_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]);
        unsafe fn mul_mat4_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]);
        unsafe fn div_mat2(&self, a: &mut [f32; 4], b: f32);
        unsafe fn div_mat3(&self, a: &mut [f32; 9], b: f32);
        unsafe fn div_mat4(&self, a: &mut [f32; 16], b: f32);
        unsafe fn div_fast_mat2(&self, a: &mut [f32; 4], b: f32);
        unsafe fn div_fast_mat3(&self, a: &mut [f32; 9], b: f32);
        unsafe fn div_fast_mat4(&self, a: &mut [f32; 16], b: f32);
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
                    b[row] = a[start_idx..end_idx].iter().zip(src_vec.iter()).map(|(x, y)| {
                        x * y
                    }).sum();
                    b[row + 2] = b[row];
                }
            }

            unsafe fn mul_vec3_mat3(&self, a: &[f32; 9], b: &mut [f32; 4]) {
                let src_vec = b.clone();
                for row in 0..3 {
                    let start_idx = row * 3;
                    let end_idx = start_idx + 3;
                    b[row] = a[start_idx..end_idx].iter().zip(src_vec.iter()).map(|(x, y)| {
                        x * y
                    }).sum();
                }
            }

            unsafe fn mul_vec4_mat4(&self, a: &[f32; 16], b: &mut [f32; 4]) {
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

            unsafe fn div_mat2(&self, a: &mut [f32; 4], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x /= b;
                });
            }

            unsafe fn div_mat3(&self, a: &mut [f32; 9], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x /= b;
                });
            }

            unsafe fn div_mat4(&self, a: &mut [f32; 16], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x /= b;
                });
            }

            unsafe fn div_fast_mat2(&self, a: &mut [f32; 4], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x /= b;
                });
            }

            unsafe fn div_fast_mat3(&self, a: &mut [f32; 9], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x /= b;
                });
            }

            unsafe fn div_fast_mat4(&self, a: &mut [f32; 16], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x /= b;
                });
            }

            unsafe fn transpose_mat2(&self, a: &mut [f32; 4]) {
                a.swap(1, 2);
            }

            unsafe fn transpose_mat3(&self, a: &mut [f32; 9]) {
                a.swap(1, 3);
                a.swap(2, 6);
                a.swap(5, 6);
            }

            unsafe fn transpose_mat4(&self, a: &mut [f32; 16]) {
                a.swap(1, 4);
                a.swap(2, 8);
                a.swap(3, 12);
                a.swap(6, 9);
                a.swap(7, 13);
                a.swap(11, 14);
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

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    pub mod x86_64 {
        #[cfg(target_arch = "x86")]
        use core::arch::x86::*;
        #[cfg(target_arch = "x86_64")]
        use core::arch::x86_64::*;

        pub const AVX_IMPL_INS: AVXImpl = AVXImpl();

        pub struct AVXImpl();

        impl super::MatImpl for AVXImpl {
            #[target_feature(enable = "sse")]
            unsafe fn add_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                let xmm_0 = [
                    _mm_load_ps(&a[0]),
                    _mm_load_ps(&b[0]),
                ];
                let xmm_1 = _mm_add_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ps(&mut a[0], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn add_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_load_ps(&b[0])
                ];
                let ymm_1 = _mm256_add_ps(ymm_0[0], ymm_0[1]);
                _mm256_store_ps(&mut a[0], ymm_1);
                let xmm_0 = [
                    _mm_load_ss(&a[8]),
                    _mm_load_ss(&b[8]),
                ];
                let xmm_1 = _mm_add_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ss(&mut a[8], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn add_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_load_ps(&a[8]),
                    _mm256_load_ps(&b[0]),
                    _mm256_load_ps(&b[8]),
                ];
                let ymm_1 = [
                    _mm256_add_ps(ymm_0[0], ymm_0[2]),
                    _mm256_add_ps(ymm_0[1], ymm_0[3]),
                ];
                _mm256_store_ps(&mut a[0], ymm_1[0]);
                _mm256_store_ps(&mut a[8], ymm_1[1]);
            }

            #[target_feature(enable = "sse")]
            unsafe fn sub_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                let xmm_0 = [
                    _mm_load_ps(&a[0]),
                    _mm_load_ps(&b[0]),
                ];
                let xmm_1 = _mm_sub_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ps(&mut a[0], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn sub_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_load_ps(&b[0])
                ];
                let ymm_1 = _mm256_sub_ps(ymm_0[0], ymm_0[1]);
                _mm256_store_ps(&mut a[0], ymm_1);
                let xmm_0 = [
                    _mm_load_ss(&a[8]),
                    _mm_load_ss(&b[8]),
                ];
                let xmm_1 = _mm_sub_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ss(&mut a[8], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn sub_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_load_ps(&a[8]),
                    _mm256_load_ps(&b[0]),
                    _mm256_load_ps(&b[8]),
                ];
                let ymm_1 = [
                    _mm256_sub_ps(ymm_0[0], ymm_0[2]),
                    _mm256_sub_ps(ymm_0[1], ymm_0[3]),
                ];
                _mm256_store_ps(&mut a[0], ymm_1[0]);
                _mm256_store_ps(&mut a[8], ymm_1[1]);
            }

            #[target_feature(enable = "sse")]
            unsafe fn mul_f32_mat2(&self, a: &mut [f32; 4], b: f32) {
                let xmm_0 = [
                    _mm_load_ps(&a[0]),
                    _mm_set1_ps(b),
                ];
                let xmm_1 = _mm_mul_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ps(&mut a[0], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn mul_f32_mat3(&self, a: &mut [f32; 9], b: f32) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_set1_ps(b),
                ];
                let ymm_1 = _mm256_mul_ps(ymm_0[0], ymm_0[1]);
                _mm256_store_ps(&mut a[0], ymm_1);
                let xmm_0 = [
                    _mm_load_ss(&a[8]),
                    _mm_load_ss(&b)
                ];
                let xmm_1 = _mm_mul_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ss(&mut a[8], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn mul_f32_mat4(&self, a: &mut [f32; 16], b: f32) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_load_ps(&a[8]),
                    _mm256_set1_ps(b),
                ];
                let ymm_1 = [
                    _mm256_mul_ps(ymm_0[0], ymm_0[2]),
                    _mm256_mul_ps(ymm_0[1], ymm_0[2]),
                ];
                _mm256_store_ps(&mut a[0], ymm_1[0]);
                _mm256_store_ps(&mut a[8], ymm_1[1]);
            }

            #[target_feature(enable = "sse")]
            unsafe fn mul_vec2_mat2(&self, a: &[f32; 4], b: &mut [f32; 4]) {
                let xmm_0 = [
                    _mm_load_ps(&a[0]),
                    _mm_load_ps(&b[0]),
                ];
                let xmm_1 = _mm_mul_ps(xmm_0[0], xmm_0[1]);
                let xmm_2 = [
                    _mm_shuffle_ps(xmm_1, xmm_1, 0x88),
                    _mm_shuffle_ps(xmm_1, xmm_1, 0xDD),
                ];
                let xmm_3 = _mm_add_ps(xmm_2[0], xmm_2[1]);
                _mm_store_ps(&mut b[0], xmm_3);
            }

            unsafe fn mul_vec3_mat3(&self, a: &[f32; 9], b: &mut [f32; 4]) {
                unimplemented!()
            }

            #[target_feature(enable = "avx")]
            unsafe fn mul_vec4_mat4(&self, a: &[f32; 16], b: &mut [f32; 4]) {
                let xmm_0 = _mm_load_ps(&b[0]);
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_load_ps(&a[8]),
                    _mm256_setr_m128(xmm_0, xmm_0),
                ];
                let ymm_1 = [
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[0], 0x08),
                        _mm256_permute_ps(ymm_0[1], 0x80),
                        0x33,
                    ),
                    _mm256_permute_ps(ymm_0[2], 0x88),
                ];
                let ymm_2 = _mm256_mul_ps(ymm_1[0], ymm_1[1]);
                let ymm_3 = [
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[0], 0x0D),
                        _mm256_permute_ps(ymm_0[1], 0xD0),
                        0x33,
                    ),
                    _mm256_permute_ps(ymm_0[2], 0xDD),
                ];
                let ymm_4 = _mm256_add_ps(_mm256_mul_ps(ymm_3[0], ymm_3[1]), ymm_2);
                let xmm_1 = _mm_add_ps(
                    _mm256_extractf128_ps(ymm_4, 0x00),
                    _mm256_extractf128_ps(ymm_4, 0x01),
                );
                _mm_store_ps(&mut b[0], xmm_1);
            }

            #[target_feature(enable = "sse")]
            unsafe fn mul_mat2_mat2(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                let xmm_0 = [
                    _mm_load_ps(&a[0]),
                    _mm_load_ps(&b[0]),
                ];
                let xmm_1 = [
                    _mm_shuffle_ps(xmm_0[0], xmm_0[0], 0xA0),
                    _mm_movelh_ps(xmm_0[1], xmm_0[1]),
                    _mm_shuffle_ps(xmm_0[0], xmm_0[0], 0xF5),
                    _mm_movehl_ps(xmm_0[1], xmm_0[1]),
                ];
                let xmm_2 = [
                    _mm_mul_ps(xmm_1[0], xmm_1[1]),
                    _mm_mul_ps(xmm_1[2], xmm_1[3]),
                ];
                let xmm_3 = _mm_add_ps(xmm_2[0], xmm_2[1]);
                _mm_store_ps(&mut a[0], xmm_3);
            }

            unsafe fn mul_mat3_mat3(&self, a: &mut [f32; 9], b: &[f32; 9]) {
                unimplemented!()
            }

            #[target_feature(enable = "avx")]
            unsafe fn mul_mat4_mat4(&self, a: &mut [f32; 16], b: &[f32; 16]) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_load_ps(&a[8]),
                    _mm256_load_ps(&b[0]),
                    _mm256_load_ps(&b[8]),
                ];
                let ymm_1 = [
                    _mm256_shuffle_ps(ymm_0[0], ymm_0[0], 0x00),
                    _mm256_shuffle_ps(ymm_0[1], ymm_0[1], 0x00),
                    _mm256_permute2f128_ps(ymm_0[2], ymm_0[2], 0x00),
                ];
                let ymm_2 = [
                    _mm256_mul_ps(ymm_1[0], ymm_1[2]),
                    _mm256_mul_ps(ymm_1[1], ymm_1[2]),
                ];
                let ymm_3 = [
                    _mm256_shuffle_ps(ymm_0[0], ymm_0[0], 0x55),
                    _mm256_shuffle_ps(ymm_0[1], ymm_0[1], 0x55),
                    _mm256_permute2f128_ps(ymm_0[2], ymm_0[2], 0x11),
                ];
                let ymm_4 = [
                    _mm256_add_ps(_mm256_mul_ps(ymm_3[0], ymm_3[2]), ymm_2[0]),
                    _mm256_add_ps(_mm256_mul_ps(ymm_3[1], ymm_3[2]), ymm_2[1]),
                ];
                let ymm_5 = [
                    _mm256_shuffle_ps(ymm_0[0], ymm_0[0], 0xAA),
                    _mm256_shuffle_ps(ymm_0[1], ymm_0[1], 0xAA),
                    _mm256_permute2f128_ps(ymm_0[3], ymm_0[3], 0x00),
                ];
                let ymm_6 = [
                    _mm256_add_ps(_mm256_mul_ps(ymm_5[0], ymm_5[2]), ymm_4[0]),
                    _mm256_add_ps(_mm256_mul_ps(ymm_5[1], ymm_5[2]), ymm_4[1]),
                ];
                let ymm_7 = [
                    _mm256_shuffle_ps(ymm_0[0], ymm_0[0], 0xFF),
                    _mm256_shuffle_ps(ymm_0[1], ymm_0[1], 0xFF),
                    _mm256_permute2f128_ps(ymm_0[3], ymm_0[3], 0x11),
                ];
                let ymm_8 = [
                    _mm256_add_ps(_mm256_mul_ps(ymm_7[0], ymm_7[2]), ymm_6[0]),
                    _mm256_add_ps(_mm256_mul_ps(ymm_7[1], ymm_7[2]), ymm_6[1]),
                ];
                _mm256_store_ps(&mut a[0], ymm_8[0]);
                _mm256_store_ps(&mut a[8], ymm_8[1]);
            }

            #[target_feature(enable = "sse")]
            unsafe fn div_mat2(&self, a: &mut [f32; 4], b: f32) {
                let xmm_0 = [
                    _mm_load_ps(&a[0]),
                    _mm_set1_ps(b),
                ];
                let xmm_1 = _mm_div_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ps(&mut a[0], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn div_mat3(&self, a: &mut [f32; 9], b: f32) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_set1_ps(b),
                ];
                let ymm_1 = _mm256_div_ps(ymm_0[0], ymm_0[1]);
                _mm256_store_ps(&mut a[0], ymm_1);
                let xmm_0 = [
                    _mm_load_ss(&a[8]),
                    _mm_load_ss(&b),
                ];
                let xmm_1 = _mm_div_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ss(&mut a[8], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn div_mat4(&self, a: &mut [f32; 16], b: f32) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_load_ps(&a[8]),
                    _mm256_set1_ps(b),
                ];
                let ymm_1 = [
                    _mm256_div_ps(ymm_0[0], ymm_0[2]),
                    _mm256_div_ps(ymm_0[1], ymm_0[2]),
                ];
                _mm256_store_ps(&mut a[0], ymm_1[0]);
                _mm256_store_ps(&mut a[8], ymm_1[1]);
            }

            #[target_feature(enable = "sse")]
            unsafe fn div_fast_mat2(&self, a: &mut [f32; 4], b: f32) {
                let xmm_0 = [
                    _mm_load_ps(&a[0]),
                    _mm_rcp_ps(_mm_set1_ps(b)),
                ];
                let xmm_1 = _mm_mul_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ps(&mut a[0], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn div_fast_mat3(&self, a: &mut [f32; 9], b: f32) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_rcp_ps(_mm256_set1_ps(b)),
                ];
                let ymm_1 = _mm256_mul_ps(ymm_0[0], ymm_0[1]);
                _mm256_store_ps(&mut a[0], ymm_1);
                let xmm_0 = [
                    _mm_load_ss(&a[8]),
                    _mm_rcp_ps(_mm_load_ss(&b)),
                ];
                let xmm_1 = _mm_mul_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ss(&mut a[8], xmm_1);
            }

            #[target_feature(enable = "avx")]
            unsafe fn div_fast_mat4(&self, a: &mut [f32; 16], b: f32) {
                let ymm_0 = [
                    _mm256_load_ps(&a[0]),
                    _mm256_load_ps(&a[8]),
                    _mm256_rcp_ps(_mm256_set1_ps(b)),
                ];
                let ymm_1 = [
                    _mm256_mul_ps(ymm_0[0], ymm_0[2]),
                    _mm256_mul_ps(ymm_0[1], ymm_0[2]),
                ];
                _mm256_store_ps(&mut a[0], ymm_1[0]);
                _mm256_store_ps(&mut a[8], ymm_1[1]);
            }

            unsafe fn transpose_mat2(&self, a: &mut [f32; 4]) {
                a.swap(1, 2);
            }

            unsafe fn transpose_mat3(&self, a: &mut [f32; 9]) {
                a.swap(1, 3);
                a.swap(2, 6);
                a.swap(5, 7);
            }

            #[target_feature(enable = "sse")]
            unsafe fn transpose_mat4(&self, a: &mut [f32; 16]) {
                let mut xmm_0 = _mm_load_ps(&a[0]);
                let mut xmm_1 = _mm_load_ps(&a[4]);
                let mut xmm_2 = _mm_load_ps(&a[8]);
                let mut xmm_3 = _mm_load_ps(&a[12]);
                _MM_TRANSPOSE4_PS(&mut xmm_0, &mut xmm_1, &mut xmm_2, &mut xmm_3);
                _mm_store_ps(&mut a[0], xmm_0);
                _mm_store_ps(&mut a[4], xmm_1);
                _mm_store_ps(&mut a[8], xmm_2);
                _mm_store_ps(&mut a[12], xmm_3);
            }

            unsafe fn det_mat2(&self, a: &[f32; 4]) -> f32 {
                unimplemented!()
            }

            unsafe fn det_mat3(&self, a: &[f32; 9]) -> f32 {
                unimplemented!()
            }

            #[target_feature(enable = "avx")]
            unsafe fn det_mat4(&self, a: &[f32; 16]) -> f32 {
                #[repr(align(32))]
                struct VecRs([f32; 8]);
                let xmm_0 = [
                    _mm_load_ps(&a[0]),
                    _mm_load_ps(&a[4]),
                    _mm_load_ps(&a[8]),
                    _mm_load_ps(&a[12]),
                ];
                let ymm_0 = [
                    _mm256_set_m128(xmm_0[0], xmm_0[0]),
                    _mm256_set_m128(xmm_0[1], xmm_0[1]),
                    _mm256_set_m128(xmm_0[2], xmm_0[2]),
                    _mm256_set_m128(xmm_0[3], xmm_0[3]),
                ];
                let ymm_1 = [
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[0], 0x40),
                        _mm256_permute_ps(ymm_0[0], 0x09),
                        0xF0,
                    ),
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[1], 0xB9),
                        _mm256_permute_ps(ymm_0[1], 0x0F),
                        0xF0,
                    ),
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[0], 0xB9),
                        _mm256_permute_ps(ymm_0[0], 0x0F),
                        0xF0,
                    ),
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[1], 0x40),
                        _mm256_permute_ps(ymm_0[1], 0x09),
                        0xF0,
                    ),
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[2], 0x16),
                        _mm256_permute_ps(ymm_0[2], 0x00),
                        0xF0,
                    ),
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[3], 0xEF),
                        _mm256_permute_ps(ymm_0[3], 0x06),
                        0xF0,
                    ),
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[2], 0xEF),
                        _mm256_permute_ps(ymm_0[2], 0x06),
                        0xF0,
                    ),
                    _mm256_blend_ps(
                        _mm256_permute_ps(ymm_0[3], 0x16),
                        _mm256_permute_ps(ymm_0[3], 0x00),
                        0xF0,
                    ),
                ];
                let ymm_2 = [
                    _mm256_mul_ps(ymm_1[0], ymm_1[1]),
                    _mm256_mul_ps(ymm_1[2], ymm_1[3]),
                    _mm256_mul_ps(ymm_1[4], ymm_1[5]),
                    _mm256_mul_ps(ymm_1[6], ymm_1[7]),
                ];
                let ymm_3 = [
                    _mm256_sub_ps(ymm_2[0], ymm_2[1]),
                    _mm256_sub_ps(ymm_2[2], ymm_2[3]),
                ];
                let ymm_4 = _mm256_mul_ps(ymm_3[0], ymm_3[1]);
                let mut vec_rs = VecRs([0.0; 8]);
                _mm256_store_ps(&mut vec_rs.0[0], ymm_4);
                vec_rs.0[0] - vec_rs.0[1] + vec_rs.0[2] + vec_rs.0[3] - vec_rs.0[4]
                    + vec_rs.0[5]
            }
        }
    }
}