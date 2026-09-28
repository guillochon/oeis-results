#!/bin/bash
# Compile the kernel to a cubin (sm_86 SASS) for the native Windows driver; the Rust program embeds
# it (rust/src/gpu.rs). A cubin needs no JIT, so a newer toolkit than the driver is fine.
cd "$(dirname "$0")"
/usr/local/cuda-13.3/bin/nvcc -O3 -arch=sm_86 -w -cubin -o a8.cubin a8.cu && echo "wrote a8.cubin ($(wc -c < a8.cubin) bytes)"
