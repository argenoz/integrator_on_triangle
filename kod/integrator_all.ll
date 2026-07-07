; ModuleID = 'integrator_all.c'
source_filename = "integrator_all.c"
target datalayout = "e-m:e-p:32:32-p10:8:8-p20:8:8-i64:64-n32:64-S128-ni:1:10:20"
target triple = "wasm32"

@alpha_beta = hidden global [3 x [2 x float]] [[2 x float] [float 0.000000e+00, float 5.000000e-01], [2 x float] [float 5.000000e-01, float 0.000000e+00], [2 x float] [float 5.000000e-01, float 5.000000e-01]], align 16
@J = hidden global float 0.000000e+00, align 4
@ax = hidden global float 0.000000e+00, align 4
@ay = hidden global float 0.000000e+00, align 4
@bx = hidden global float 0.000000e+00, align 4
@by = hidden global float 0.000000e+00, align 4
@cx = hidden global float 0.000000e+00, align 4
@cy = hidden global float 0.000000e+00, align 4
@acx = hidden global float 0.000000e+00, align 4
@bcx = hidden global float 0.000000e+00, align 4
@acy = hidden global float 0.000000e+00, align 4
@bcy = hidden global float 0.000000e+00, align 4
@optimuz = hidden global [16 x float] [float 0x3FBE4BCD80000000, float 0x3FBE4BCD80000000, float 0xBFC0D4C780000000, float 0x3FD7904A80000000, float 0xBF8FBB7020000000, float 0xBF8FBB7020000000, float 0x3F9F82B260000000, float 0x3F92E9E8C0000000, float 0x3F92E9E8C0000000, float 0x3F9F8C2740000000, float 0x3FB4854860000000, float 0x3F9F8C2740000000, float 0x3FB4854860000000, float 0xBFC688AE60000000, float 0xBFC688AE60000000, float 0x3FBEBC1A40000000], align 16
@ab_pairs = hidden global [16 x [2 x float]] [[2 x float] [float 0x3FE5555560000000, float 0x3FC5555560000000], [2 x float] [float 0x3FC5555560000000, float 0x3FE5555560000000], [2 x float] [float 0x3FC5555560000000, float 0x3FC5555560000000], [2 x float] [float 0x3FD5555560000000, float 0x3FD5555560000000], [2 x float] [float 5.000000e-01, float 0.000000e+00], [2 x float] [float 0.000000e+00, float 5.000000e-01], [2 x float] [float 5.000000e-01, float 5.000000e-01], [2 x float] [float 2.500000e-01, float 7.500000e-01], [2 x float] [float 7.500000e-01, float 2.500000e-01], [2 x float] [float 0.000000e+00, float 7.500000e-01], [2 x float] [float 0.000000e+00, float 2.500000e-01], [2 x float] [float 7.500000e-01, float 0.000000e+00], [2 x float] [float 2.500000e-01, float 0.000000e+00], [2 x float] [float 2.500000e-01, float 5.000000e-01], [2 x float] [float 5.000000e-01, float 2.500000e-01], [2 x float] [float 2.500000e-01, float 2.500000e-01]], align 16

; Function Attrs: noinline nounwind optnone
define hidden float @integrate_p3() #0 {
  %1 = load float, float* @J, align 4
  %2 = load float, float* @ax, align 4
  %3 = load float, float* @ay, align 4
  %4 = call float @f(float noundef %2, float noundef %3)
  %5 = load float, float* @bx, align 4
  %6 = load float, float* @by, align 4
  %7 = call float @f(float noundef %5, float noundef %6)
  %8 = fadd float %4, %7
  %9 = load float, float* @cx, align 4
  %10 = load float, float* @cy, align 4
  %11 = call float @f(float noundef %9, float noundef %10)
  %12 = fadd float %8, %11
  %13 = fmul float %1, %12
  %14 = fdiv float %13, 6.000000e+00
  ret float %14
}

declare float @f(float noundef, float noundef) #1

; Function Attrs: noinline nounwind optnone
define hidden float @integrate_p6() #0 {
  %1 = alloca float, align 4
  %2 = alloca float, align 4
  %3 = alloca float, align 4
  %4 = load float, float* @acx, align 4
  %5 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 0, i32 0), align 16
  %6 = load float, float* @bcx, align 4
  %7 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 0, i32 1), align 4
  %8 = fmul float %6, %7
  %9 = call float @llvm.fmuladd.f32(float %4, float %5, float %8)
  %10 = load float, float* @cx, align 4
  %11 = fadd float %9, %10
  %12 = load float, float* @acy, align 4
  %13 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 0, i32 0), align 16
  %14 = load float, float* @bcy, align 4
  %15 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 0, i32 1), align 4
  %16 = fmul float %14, %15
  %17 = call float @llvm.fmuladd.f32(float %12, float %13, float %16)
  %18 = load float, float* @cy, align 4
  %19 = fadd float %17, %18
  %20 = call float @f(float noundef %11, float noundef %19)
  store float %20, float* %1, align 4
  %21 = load float, float* @acx, align 4
  %22 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 1, i32 0), align 8
  %23 = load float, float* @bcx, align 4
  %24 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 1, i32 1), align 4
  %25 = fmul float %23, %24
  %26 = call float @llvm.fmuladd.f32(float %21, float %22, float %25)
  %27 = load float, float* @cx, align 4
  %28 = fadd float %26, %27
  %29 = load float, float* @acy, align 4
  %30 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 1, i32 0), align 8
  %31 = load float, float* @bcy, align 4
  %32 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 1, i32 1), align 4
  %33 = fmul float %31, %32
  %34 = call float @llvm.fmuladd.f32(float %29, float %30, float %33)
  %35 = load float, float* @cy, align 4
  %36 = fadd float %34, %35
  %37 = call float @f(float noundef %28, float noundef %36)
  store float %37, float* %2, align 4
  %38 = load float, float* @acx, align 4
  %39 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 2, i32 0), align 16
  %40 = load float, float* @bcx, align 4
  %41 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 2, i32 1), align 4
  %42 = fmul float %40, %41
  %43 = call float @llvm.fmuladd.f32(float %38, float %39, float %42)
  %44 = load float, float* @cx, align 4
  %45 = fadd float %43, %44
  %46 = load float, float* @acy, align 4
  %47 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 2, i32 0), align 16
  %48 = load float, float* @bcy, align 4
  %49 = load float, float* getelementptr inbounds ([3 x [2 x float]], [3 x [2 x float]]* @alpha_beta, i32 0, i32 2, i32 1), align 4
  %50 = fmul float %48, %49
  %51 = call float @llvm.fmuladd.f32(float %46, float %47, float %50)
  %52 = load float, float* @cy, align 4
  %53 = fadd float %51, %52
  %54 = call float @f(float noundef %45, float noundef %53)
  store float %54, float* %3, align 4
  %55 = load float, float* %1, align 4
  %56 = load float, float* %2, align 4
  %57 = fadd float %55, %56
  %58 = load float, float* %3, align 4
  %59 = fadd float %57, %58
  %60 = load float, float* @J, align 4
  %61 = fdiv float %60, 6.000000e+00
  %62 = fmul float %59, %61
  ret float %62
}

; Function Attrs: nofree nosync nounwind readnone speculatable willreturn
declare float @llvm.fmuladd.f32(float, float, float) #2

; Function Attrs: noinline nounwind optnone
define hidden void @setV(float noundef %0, float noundef %1, float noundef %2, float noundef %3, float noundef %4, float noundef %5) #0 {
  %7 = alloca float, align 4
  %8 = alloca float, align 4
  %9 = alloca float, align 4
  %10 = alloca float, align 4
  %11 = alloca float, align 4
  %12 = alloca float, align 4
  store float %0, float* %7, align 4
  store float %1, float* %8, align 4
  store float %2, float* %9, align 4
  store float %3, float* %10, align 4
  store float %4, float* %11, align 4
  store float %5, float* %12, align 4
  %13 = load float, float* %7, align 4
  store float %13, float* @ax, align 4
  %14 = load float, float* %8, align 4
  store float %14, float* @ay, align 4
  %15 = load float, float* %9, align 4
  store float %15, float* @bx, align 4
  %16 = load float, float* %10, align 4
  store float %16, float* @by, align 4
  %17 = load float, float* %11, align 4
  store float %17, float* @cx, align 4
  %18 = load float, float* %12, align 4
  store float %18, float* @cy, align 4
  %19 = load float, float* @ax, align 4
  %20 = load float, float* @cx, align 4
  %21 = fsub float %19, %20
  store float %21, float* @acx, align 4
  %22 = load float, float* @ay, align 4
  %23 = load float, float* @cy, align 4
  %24 = fsub float %22, %23
  store float %24, float* @acy, align 4
  %25 = load float, float* @bx, align 4
  %26 = load float, float* @cx, align 4
  %27 = fsub float %25, %26
  store float %27, float* @bcx, align 4
  %28 = load float, float* @by, align 4
  %29 = load float, float* @cy, align 4
  %30 = fsub float %28, %29
  store float %30, float* @bcy, align 4
  %31 = load float, float* @acx, align 4
  %32 = load float, float* @bcy, align 4
  %33 = load float, float* @acy, align 4
  %34 = load float, float* @bcx, align 4
  %35 = fmul float %33, %34
  %36 = fneg float %35
  %37 = call float @llvm.fmuladd.f32(float %31, float %32, float %36)
  store float %37, float* @J, align 4
  ret void
}

; Function Attrs: noinline nounwind optnone
define hidden float @integrate_p16() #0 {
  %1 = alloca i32, align 4
  %2 = alloca i32, align 4
  %3 = alloca float, align 4
  %4 = alloca float, align 4
  %5 = alloca float, align 4
  %6 = alloca float, align 4
  %7 = alloca float, align 4
  store i32 15, i32* %1, align 4
  %8 = load i32, i32* %1, align 4
  %9 = getelementptr inbounds [16 x [2 x float]], [16 x [2 x float]]* @ab_pairs, i32 0, i32 %8
  %10 = getelementptr inbounds [2 x float], [2 x float]* %9, i32 0, i32 0
  %11 = load float, float* %10, align 8
  store float %11, float* %4, align 4
  %12 = load i32, i32* %1, align 4
  %13 = getelementptr inbounds [16 x [2 x float]], [16 x [2 x float]]* @ab_pairs, i32 0, i32 %12
  %14 = getelementptr inbounds [2 x float], [2 x float]* %13, i32 0, i32 1
  %15 = load float, float* %14, align 4
  store float %15, float* %5, align 4
  %16 = load float, float* @acx, align 4
  %17 = load float, float* %4, align 4
  %18 = load float, float* @bcx, align 4
  %19 = load float, float* %5, align 4
  %20 = fmul float %18, %19
  %21 = call float @llvm.fmuladd.f32(float %16, float %17, float %20)
  %22 = load float, float* @cx, align 4
  %23 = fadd float %21, %22
  store float %23, float* %6, align 4
  %24 = load float, float* @acy, align 4
  %25 = load float, float* %4, align 4
  %26 = load float, float* @bcy, align 4
  %27 = load float, float* %5, align 4
  %28 = fmul float %26, %27
  %29 = call float @llvm.fmuladd.f32(float %24, float %25, float %28)
  %30 = load float, float* @cy, align 4
  %31 = fadd float %29, %30
  store float %31, float* %7, align 4
  %32 = load i32, i32* %1, align 4
  %33 = getelementptr inbounds [16 x float], [16 x float]* @optimuz, i32 0, i32 %32
  %34 = load float, float* %33, align 4
  %35 = load float, float* %6, align 4
  %36 = load float, float* %7, align 4
  %37 = call float @f(float noundef %35, float noundef %36)
  %38 = fmul float %34, %37
  store float %38, float* %3, align 4
  store i32 15, i32* %1, align 4
  br label %39

39:                                               ; preds = %74, %0
  %40 = load i32, i32* %1, align 4
  %41 = add nsw i32 %40, -1
  store i32 %41, i32* %1, align 4
  %42 = load i32, i32* %1, align 4
  %43 = getelementptr inbounds [16 x [2 x float]], [16 x [2 x float]]* @ab_pairs, i32 0, i32 %42
  %44 = getelementptr inbounds [2 x float], [2 x float]* %43, i32 0, i32 0
  %45 = load float, float* %44, align 8
  store float %45, float* %4, align 4
  %46 = load i32, i32* %1, align 4
  %47 = getelementptr inbounds [16 x [2 x float]], [16 x [2 x float]]* @ab_pairs, i32 0, i32 %46
  %48 = getelementptr inbounds [2 x float], [2 x float]* %47, i32 0, i32 1
  %49 = load float, float* %48, align 4
  store float %49, float* %5, align 4
  %50 = load float, float* @acx, align 4
  %51 = load float, float* %4, align 4
  %52 = load float, float* @bcx, align 4
  %53 = load float, float* %5, align 4
  %54 = fmul float %52, %53
  %55 = call float @llvm.fmuladd.f32(float %50, float %51, float %54)
  %56 = load float, float* @cx, align 4
  %57 = fadd float %55, %56
  store float %57, float* %6, align 4
  %58 = load float, float* @acy, align 4
  %59 = load float, float* %4, align 4
  %60 = load float, float* @bcy, align 4
  %61 = load float, float* %5, align 4
  %62 = fmul float %60, %61
  %63 = call float @llvm.fmuladd.f32(float %58, float %59, float %62)
  %64 = load float, float* @cy, align 4
  %65 = fadd float %63, %64
  store float %65, float* %7, align 4
  %66 = load i32, i32* %1, align 4
  %67 = getelementptr inbounds [16 x float], [16 x float]* @optimuz, i32 0, i32 %66
  %68 = load float, float* %67, align 4
  %69 = load float, float* %6, align 4
  %70 = load float, float* %7, align 4
  %71 = call float @f(float noundef %69, float noundef %70)
  %72 = load float, float* %3, align 4
  %73 = call float @llvm.fmuladd.f32(float %68, float %71, float %72)
  store float %73, float* %3, align 4
  br label %74

74:                                               ; preds = %39
  %75 = load i32, i32* %1, align 4
  %76 = icmp ne i32 %75, 0
  br i1 %76, label %39, label %77, !llvm.loop !2

77:                                               ; preds = %74
  %78 = load float, float* %3, align 4
  %79 = load float, float* @J, align 4
  %80 = fmul float %78, %79
  ret float %80
}

attributes #0 = { noinline nounwind optnone "frame-pointer"="none" "min-legal-vector-width"="0" "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="generic" }
attributes #1 = { "frame-pointer"="none" "no-trapping-math"="true" "stack-protector-buffer-size"="8" "target-cpu"="generic" }
attributes #2 = { nofree nosync nounwind readnone speculatable willreturn }

!llvm.module.flags = !{!0}
!llvm.ident = !{!1}

!0 = !{i32 1, !"wchar_size", i32 4}
!1 = !{!"Ubuntu clang version 14.0.0-1ubuntu1.1"}
!2 = distinct !{!2, !3}
!3 = !{!"llvm.loop.mustprogress"}
