// Compiled to SPIR-V with:
// `dxc -T cs_6_0 -E "add" -spirv -fvk-use-dx-layout -fspv-target-env=vulkan1.1 -Fo "add.spv" "add.hlsl"`

[[vk::binding(0, 0)]] RWStructuredBuffer<int> OutBuffer;
[[vk::binding(1, 0)]] RWStructuredBuffer<int> InBuffer1;
[[vk::binding(2, 0)]] RWStructuredBuffer<int> InBuffer2;

[numthreads(1, 1, 1)]
void add(uint3 DTid : SV_DispatchThreadID) {
  OutBuffer[DTid.x] = InBuffer1[DTid.x] + InBuffer2[DTid.x];
}
