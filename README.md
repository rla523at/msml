# 1. Build msml
## 1.1. CMake 설치
[홈페이지](https://cmake.org/download/)에 들어가서 platform에 맞는 Binary distribution을 다운받고 실행시키면 설치가 완료된다.

## 1.2. CMake 실행
1. CMake(cmake-gui)를 실행
2. Where is the source code에 msml/library build 폴더 선택
3. msml/library build/bin 폴더 선택
4. Configure 선택 >> Specify the generator for this project >> library를 사용하고자 하는 플랫폼 선택 >> Finish
5. Generate

## 1.3 Build library
### Visual studio의 경우
1. library build >> bin >> msml.sln 실행
2. debug x64환경에서 빌드
3. relaese x64환경에서 빌드
4. msml/library build/lib 폴더에 `msml_release.lib`, `msml_debug.lib` 생성 완료

# 2. Initial Setting
1. 사용하고자 하는 프로젝트의 `.vcxproj`가 있는 위치에 lib 폴더 생성후 msml 폴더 생성
2. 생성한 project/lib/msml 폴더안에 msml/library build/lib 폴더에 있는 `msml_release.lib`, `msml_debug.lib` 복사
3. lib/msml 폴더안에 msml/src 폴더를 복사
4. lib/msml 폴더안에 미리 받은 mkl_core.lib, mkl_intel_lp64.lib, mkl_sequential.lib 복사

project/lib 폴더관련 트리구조는 다음과 같아야 한다.
```
└─lib
    └─msml
        │  mkl_core.lib
        │  mkl_intel_lp64.lib
        │  mkl_sequential.lib
        │  msml_debug.lib
        │  msml_release.lib
        │
        └─src
                Matrix.h
                Real_Vector.h
```

# 3. Visual Studio Setting Guide

1. 프로젝트 >> 속성 >> 일반 >> 플랫폼 도구 집합 >> C++17 이상
2. 프로젝트 >> 속성 >> 링커 >> 일반 >> 추가 라이브러리 디렉터리 >> 아래 내용 입력
```
$(ProjectDir)lib/msml
```
3.  프로젝트 >> 속성 >> 링커 >> 입력 >> 추가종속성 >> 편집 >> 아래 내용 입력


```
(구성: debug)
mkl_core.lib
mkl_intel_lp64.lib
mkl_sequential.lib
msml_debug.lib


(구성: release)
mkl_core.lib
mkl_intel_lp64.lib
mkl_sequential.lib
msml_release.lib
```

# 사용법
사용하고자 하는 cpp 파일에서 아래 내용 입력후 사용
```
(Real Vector)
#include "lib/msml/src/Real_Vector.h"

(Matrix)
#include "lib/msml/src/Matrix.h"
```