# Third-party tools used

This software was implemented using the Julia language, which is licensed under the [MIT
License](https://github.com/JuliaLang/julia/blob/master/LICENSE.md), and depends 
itself on other software which are listed [here](https://github.com/JuliaLang/julia/blob/master/THIRDPARTY.md).
The software is compiled using PackageCompiler.jl, which is shipped under the 
[MIT License](https://github.com/JuliaLang/PackageCompiler.jl/blob/master/LICENSE).

The compiled binary contains external libraries which are used by the standard library
of the Julia language. The following external libraries are shipped alongside the 
software:

 - [LIBBLASTRAMPOLINE](https://github.com/staticfloat/libblastrampoline/blob/main/LICENSE): MIT License        
 - [DSFMT][https://github.com/MersenneTwister-Lab/dSFMT/blob/master/LICENSE.txt]: BSD-3 License                 
 - [GCC](https://github.com/gcc-mirror/gcc/blob/master/COPYING.LIB): LGPL2.1 License and Notice
 - [GMP](https://gmplib.org/manual/Copying.html#Copying): LGPL3+ License or GPL2+ License    
 - [LLVM](https://releases.llvm.org/3.9.0/LICENSE.TXT): UIUC License                
 - [MPFR](https://www.mpfr.org/mpfr-current/mpfr.html#Copying): LGPL3+ License
 - [LIBUNWIND](https://github.com/libunwind/libunwind/blob/master/LICENSE): MIT License      
 - [LIBUV](https://github.com/JuliaLang/libuv/blob/julia-uv2-1.39.0/LICENSE): MIT License
 - [ZLIB](https://zlib.net/zlib_license.html): ZLib License
 - [OPENBLAS](https://raw.github.com/xianyi/OpenBLAS/master/LICENSE): BSD-3 License  
 - [OPENLIBM][https://github.com/JuliaMath/openlibm/blob/master/LICENSE.md]: MIT License and BSD-2 License and ISC License 
 - [PCRE](https://www.pcre.org/licence.txt): BSD-3 License

