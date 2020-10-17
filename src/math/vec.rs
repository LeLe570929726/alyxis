#[repr(align(16))]
#[derive(Clone, PartialEq, Debug)]
pub struct Vec2([f32; 4]);

#[repr(align(16))]
#[derive(Clone, PartialEq, Debug)]
pub struct Vec3([f32; 4]);

#[repr(align(16))]
#[derive(Clone, PartialEq, Debug)]
pub struct Vec4([f32; 4]);

impl Vec2 {
    pub fn new(x: f32, y: f32) -> Self {
        Self([x, y, x, y])
    }
}

impl Vec3 {
    pub fn new(x: f32, y: f32, z: f32) -> Self {
        Self([x, y, z, 0.0])
    }
}

impl Vec4 {
    pub fn new(x: f32, y: f32, z: f32, w: f32) -> Self {
        Self([x, y, z, w])
    }
}

mod impls {
    pub static mut IMPL_INS: &dyn VecImpl = &general::FPU_IMPL_INS;

    pub trait VecImpl {
        unsafe fn add(&self, a: &mut [f32; 4], b: &[f32; 4]);
        unsafe fn sub(&self, a: &mut [f32; 4], b: &[f32; 4]);
        unsafe fn mul(&self, a: &mut [f32; 4], b: f32);
        unsafe fn div(&self, a: &mut [f32; 4], b: f32);
        unsafe fn dot_vec2(&self, a: &[f32; 4], b: &[f32; 4]) -> f32;
        unsafe fn dot_vec3(&self, a: &[f32; 4], b: &[f32; 4]) -> f32;
        unsafe fn dot_vec4(&self, a: &[f32; 4], b: &[f32; 4]) -> f32;
        unsafe fn norm_vec2(&self, a: &[f32; 4]) -> f32;
        unsafe fn norm_vec3(&self, a: &[f32; 4]) -> f32;
        unsafe fn norm_vec4(&self, a: &[f32; 4]) -> f32;
        unsafe fn normalize_vec2(&self, a: &mut [f32; 4]);
        unsafe fn normalize_vec3(&self, a: &mut [f32; 4]);
        unsafe fn normalize_vec4(&self, a: &mut [f32; 4]);
    }

    pub mod general {
        pub const FPU_IMPL_INS: FPUImpl = FPUImpl();

        pub struct FPUImpl();

        impl super::VecImpl for FPUImpl {
            unsafe fn add(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
                    *x += y;
                });
            }

            unsafe fn sub(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                a.iter_mut().zip(b.iter()).for_each(|(x, y)| {
                    *x -= y;
                });
            }

            unsafe fn mul(&self, a: &mut [f32; 4], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x *= b;
                });
            }

            unsafe fn div(&self, a: &mut [f32; 4], b: f32) {
                a.iter_mut().for_each(|x| {
                    *x /= b;
                });
            }

            unsafe fn dot_vec2(&self, a: &[f32; 4], b: &[f32; 4]) -> f32 {
                a[..2].iter().zip(b.iter()).map(|(x, y)| {
                    x * y
                }).sum()
            }

            unsafe fn dot_vec3(&self, a: &[f32; 4], b: &[f32; 4]) -> f32 {
                a[..3].iter().zip(b.iter()).map(|(x, y)| {
                    x * y
                }).sum()
            }

            unsafe fn dot_vec4(&self, a: &[f32; 4], b: &[f32; 4]) -> f32 {
                a.iter().zip(b.iter()).map(|(x, y)| {
                    x * y
                }).sum()
            }

            unsafe fn norm_vec2(&self, a: &[f32; 4]) -> f32 {
                f32::sqrt(a[..2].iter().map(|x| {
                    x * x
                }).sum())
            }

            unsafe fn norm_vec3(&self, a: &[f32; 4]) -> f32 {
                f32::sqrt(a[..3].iter().map(|x| {
                    x * x
                }).sum())
            }

            unsafe fn norm_vec4(&self, a: &[f32; 4]) -> f32 {
                f32::sqrt(a.iter().map(|x| {
                    x * x
                }).sum())
            }

            unsafe fn normalize_vec2(&self, a: &mut [f32; 4]) {
                let norm = self.norm_vec2(a);
                a.iter_mut().for_each(|x| {
                    *x /= norm;
                });
            }

            unsafe fn normalize_vec3(&self, a: &mut [f32; 4]) {
                let norm = self.norm_vec3(a);
                a[..3].iter_mut().for_each(|x| {
                    *x /= norm;
                });
            }

            unsafe fn normalize_vec4(&self, a: &mut [f32; 4]) {
                let norm = self.norm_vec4(a);
                a.iter_mut().for_each(|x| {
                    *x /= norm;
                });
            }
        }
    }

    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    pub mod x86_64 {
        #[cfg(target_arch = "x86")]
        use core::arch::x86::*;
        #[cfg(target_arch = "x86_64")]
        use core::arch::x86_64::*;

        pub const SSE_IMPL_INS: SSEImpl = SSEImpl();

        pub struct SSEImpl();

        impl super::VecImpl for SSEImpl {
            #[target_feature(enable = "sse")]
            unsafe fn add(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                let xmm_0 = [
                    _mm_load_ps(a.as_ptr()),
                    _mm_load_ps(b.as_ptr()),
                ];
                let xmm_1 = _mm_add_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ps(a.as_mut_ptr(), xmm_1);
            }

            #[target_feature(enable = "sse")]
            unsafe fn sub(&self, a: &mut [f32; 4], b: &[f32; 4]) {
                let xmm_0 = [
                    _mm_load_ps(a.as_ptr()),
                    _mm_load_ps(b.as_ptr()),
                ];
                let xmm_1 = _mm_sub_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ps(a.as_mut_ptr(), xmm_1);
            }

            #[target_feature(enable = "sse")]
            unsafe fn mul(&self, a: &mut [f32; 4], b: f32) {
                let xmm_0 = [
                    _mm_load_ps(a.as_ptr()),
                    _mm_set1_ps(b),
                ];
                let xmm_1 = _mm_mul_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ps(a.as_mut_ptr(), xmm_1);
            }

            #[target_feature(enable = "sse")]
            unsafe fn div(&self, a: &mut [f32; 4], b: f32) {
                let xmm_0 = [
                    _mm_load_ps(a.as_ptr()),
                    _mm_set1_ps(b),
                ];
                let xmm_1 = _mm_div_ps(xmm_0[0], xmm_0[1]);
                _mm_store_ps(a.as_mut_ptr(), xmm_1);
            }

            #[target_feature(enable = "sse4.1")]
            unsafe fn dot_vec2(&self, a: &[f32; 4], b: &[f32; 4]) -> f32 {
                let xmm_0 = [
                    _mm_load_ps(a.as_ptr()),
                    _mm_load_ps(b.as_ptr()),
                ];
                let xmm_1 = _mm_dp_ps(xmm_0[0], xmm_0[1], 0x31);
                let mut rs = 0.0;
                _mm_store_ss(&mut rs, xmm_1);
                rs
            }

            #[target_feature(enable = "sse4.1")]
            unsafe fn dot_vec3(&self, a: &[f32; 4], b: &[f32; 4]) -> f32 {
                let xmm_0 = [
                    _mm_load_ps(a.as_ptr()),
                    _mm_load_ps(b.as_ptr()),
                ];
                let xmm_1 = _mm_dp_ps(xmm_0[0], xmm_0[1], 0x71);
                let mut rs = 0.0;
                _mm_store_ss(&mut rs, xmm_1);
                rs
            }

            #[target_feature(enable = "sse4.1")]
            unsafe fn dot_vec4(&self, a: &[f32; 4], b: &[f32; 4]) -> f32 {
                let xmm_0 = [
                    _mm_load_ps(a.as_ptr()),
                    _mm_load_ps(b.as_ptr()),
                ];
                let xmm_1 = _mm_dp_ps(xmm_0[0], xmm_0[1], 0xF1);
                let mut rs = 0.0;
                _mm_store_ss(&mut rs, xmm_1);
                rs
            }

            #[target_feature(enable = "sse4.1")]
            unsafe fn norm_vec2(&self, a: &[f32; 4]) -> f32 {
                let xmm_0 = _mm_load_ps(a.as_ptr());
                let xmm_1 = _mm_dp_ps(xmm_0, xmm_0, 0x33);
                let xmm_2 = _mm_sqrt_ps(xmm_1);
                let mut rs = 0.0;
                _mm_store_ss(&mut rs, xmm_2);
                rs
            }

            #[target_feature(enable = "sse4.1")]
            unsafe fn norm_vec3(&self, a: &[f32; 4]) -> f32 {
                let xmm_0 = _mm_load_ps(a.as_ptr());
                let xmm_1 = _mm_dp_ps(xmm_0, xmm_0, 0x77);
                let xmm_2 = _mm_sqrt_ps(xmm_1);
                let mut rs = 0.0;
                _mm_store_ss(&mut rs, xmm_2);
                rs
            }

            #[target_feature(enable = "sse4.1")]
            unsafe fn norm_vec4(&self, a: &[f32; 4]) -> f32 {
                let xmm_0 = _mm_load_ps(a.as_ptr());
                let xmm_1 = _mm_dp_ps(xmm_0, xmm_0, 0xFF);
                let xmm_2 = _mm_sqrt_ps(xmm_1);
                let mut rs = 0.0;
                _mm_store_ss(&mut rs, xmm_2);
                rs
            }

            #[target_feature(enable = "sse4.1")]
            unsafe fn normalize_vec2(&self, a: &mut [f32; 4]) {
                let xmm_0 = _mm_load_ps(a.as_ptr());
                let xmm_1 = _mm_dp_ps(xmm_0, xmm_0, 0x3F);
                let xmm_2 = _mm_rsqrt_ps(xmm_1);
                let xmm_3 = _mm_mul_ps(xmm_0, xmm_2);
                _mm_store_ps(a.as_mut_ptr(), xmm_3);
            }

            #[target_feature(enable = "sse4.1")]
            unsafe fn normalize_vec3(&self, a: &mut [f32; 4]) {
                let xmm_0 = _mm_load_ps(a.as_ptr());
                let xmm_1 = _mm_dp_ps(xmm_0, xmm_0, 0x77);
                let xmm_2 = _mm_rsqrt_ps(xmm_1);
                let xmm_3 = _mm_mul_ps(xmm_0, xmm_2);
                _mm_store_ps(a.as_mut_ptr(), xmm_3);
            }

            #[target_feature(enable = "sse4.1")]
            unsafe fn normalize_vec4(&self, a: &mut [f32; 4]) {
                let xmm_0 = _mm_load_ps(a.as_ptr());
                let xmm_1 = _mm_dp_ps(xmm_0, xmm_0, 0xFF);
                let xmm_2 = _mm_rsqrt_ps(xmm_1);
                let xmm_3 = _mm_mul_ps(xmm_0, xmm_2);
                _mm_store_ps(a.as_mut_ptr(), xmm_3);
            }
        }
    }
}

#[cfg(test)]
mod test {
    const VEC2_A: super::Vec2 = super::Vec2([1.0, 2.0, 1.0, 2.0]);
    const VEC2_B: super::Vec2 = super::Vec2([-3.0, -4.0, -3.0, -4.0]);
    const VEC3_A: super::Vec3 = super::Vec3([1.0, 2.0, 3.0, 0.0]);
    const VEC3_B: super::Vec3 = super::Vec3([-3.0, -2.0, -1.0, 0.0]);
    const VEC4_A: super::Vec4 = super::Vec4([1.0, 2.0, 3.0, 4.0]);
    const VEC4_B: super::Vec4 = super::Vec4([-4.0, -3.0, -2.0, -1.0]);

    mod impls {
        mod general {
            use super::super::super::*;
            use super::super::super::impls::VecImpl;

            #[test]
            fn add() {
                let mut vec_rs = test::VEC4_A.clone();
                let vec_rs_t = [-3.0, -1.0, 1.0, 3.0];
                unsafe { impls::general::FPU_IMPL_INS.add(&mut vec_rs.0, &test::VEC4_B.0); }
                assert_eq!(vec_rs.0, vec_rs_t);
            }

            #[test]
            fn sub() {
                let mut vec_rs = test::VEC4_A.clone();
                let vec_rs_t = [5.0, 5.0, 5.0, 5.0];
                unsafe { impls::general::FPU_IMPL_INS.sub(&mut vec_rs.0, &test::VEC4_B.0); }
                assert_eq!(vec_rs.0, vec_rs_t);
            }

            #[test]
            fn mul() {
                let mut vec_rs = test::VEC4_A.clone();
                let vec_rs_t = [2.0, 4.0, 6.0, 8.0];
                unsafe { impls::general::FPU_IMPL_INS.mul(&mut vec_rs.0, 2.0); }
                assert_eq!(vec_rs.0, vec_rs_t);
            }

            #[test]
            fn div() {
                let mut vec_rs = test::VEC4_A.clone();
                let vec_rs_t = [0.5, 1.0, 1.5, 2.0];
                unsafe { impls::general::FPU_IMPL_INS.div(&mut vec_rs.0, 2.0); }
                assert_eq!(vec_rs.0, vec_rs_t);
            }

            #[test]
            fn dot_vec2() {
                let f32_rs = unsafe { impls::general::FPU_IMPL_INS.dot_vec2(&test::VEC2_A.0, &test::VEC2_B.0) };
                let f32_rs_t = -11.0;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn dot_vec3() {
                let f32_rs = unsafe { impls::general::FPU_IMPL_INS.dot_vec3(&test::VEC3_A.0, &test::VEC3_B.0) };
                let f32_rs_t = -10.0;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn dot_vec4() {
                let f32_rs = unsafe { impls::general::FPU_IMPL_INS.dot_vec4(&test::VEC4_A.0, &test::VEC4_B.0) };
                let f32_rs_t = -20.0;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn norm_vec2() {
                let f32_rs = unsafe { impls::general::FPU_IMPL_INS.norm_vec2(&test::VEC2_A.0) };
                let f32_rs_t = 2.236068;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn norm_vec3() {
                let f32_rs = unsafe { impls::general::FPU_IMPL_INS.norm_vec3(&test::VEC3_A.0) };
                let f32_rs_t = 3.7416575;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn norm_vec4() {
                let f32_rs = unsafe { impls::general::FPU_IMPL_INS.norm_vec4(&test::VEC4_A.0) };
                let f32_rs_t = 5.477226;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn normalize_vec2() {
                let mut vec_rs = test::VEC2_A.clone();
                let vec_rs_t = [0.4472136, 0.8944272, 0.4472136, 0.8944272];
                unsafe { impls::general::FPU_IMPL_INS.normalize_vec2(&mut vec_rs.0) };
                assert_eq!(vec_rs.0, vec_rs_t);
            }
        }

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        mod x86_64 {
            use super::super::super::*;
            use super::super::super::impls::VecImpl;

            #[test]
            fn add() {
                let mut vec_rs = test::VEC4_A.clone();
                let vec_rs_t = [-3.0, -1.0, 1.0, 3.0];
                unsafe { impls::x86_64::SSE_IMPL_INS.add(&mut vec_rs.0, &test::VEC4_B.0); }
                assert_eq!(vec_rs.0, vec_rs_t);
            }

            #[test]
            fn sub() {
                let mut vec_rs = test::VEC4_A.clone();
                let vec_rs_t = [5.0, 5.0, 5.0, 5.0];
                unsafe { impls::x86_64::SSE_IMPL_INS.sub(&mut vec_rs.0, &test::VEC4_B.0); }
                assert_eq!(vec_rs.0, vec_rs_t);
            }

            #[test]
            fn mul() {
                let mut vec_rs = test::VEC4_A.clone();
                let vec_rs_t = [2.0, 4.0, 6.0, 8.0];
                unsafe { impls::x86_64::SSE_IMPL_INS.mul(&mut vec_rs.0, 2.0); }
                assert_eq!(vec_rs.0, vec_rs_t);
            }

            #[test]
            fn div() {
                let mut vec_rs = test::VEC4_A.clone();
                let vec_rs_t = [0.5, 1.0, 1.5, 2.0];
                unsafe { impls::x86_64::SSE_IMPL_INS.div(&mut vec_rs.0, 2.0); }
                assert_eq!(vec_rs.0, vec_rs_t);
            }

            #[test]
            fn dot_vec2() {
                let f32_rs = unsafe { impls::x86_64::SSE_IMPL_INS.dot_vec2(&test::VEC2_A.0, &test::VEC2_B.0) };
                let f32_rs_t = -11.0;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn dot_vec3() {
                let f32_rs = unsafe { impls::x86_64::SSE_IMPL_INS.dot_vec3(&test::VEC3_A.0, &test::VEC3_B.0) };
                let f32_rs_t = -10.0;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn dot_vec4() {
                let f32_rs = unsafe { impls::x86_64::SSE_IMPL_INS.dot_vec4(&test::VEC4_A.0, &test::VEC4_B.0) };
                let f32_rs_t = -20.0;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn norm_vec2() {
                let f32_rs = unsafe { impls::x86_64::SSE_IMPL_INS.norm_vec2(&test::VEC2_A.0) };
                let f32_rs_t = 2.236068;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn norm_vec3() {
                let f32_rs = unsafe { impls::x86_64::SSE_IMPL_INS.norm_vec3(&test::VEC3_A.0) };
                let f32_rs_t = 3.7416575;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn norm_vec4() {
                let f32_rs = unsafe { impls::x86_64::SSE_IMPL_INS.norm_vec4(&test::VEC4_A.0) };
                let f32_rs_t = 5.477226;
                assert_eq!(f32_rs, f32_rs_t);
            }

            #[test]
            fn normalize_vec2() {
                let mut vec_rs = test::VEC2_A.clone();
                let vec_rs_t = [0.44714355, 0.8942871, 0.44714355, 0.8942871];
                unsafe { impls::x86_64::SSE_IMPL_INS.normalize_vec2(&mut vec_rs.0) };
                assert_eq!(vec_rs.0, vec_rs_t);
            }
        }
    }
}