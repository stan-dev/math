# run from the build folder
import asyncio
import boost_cobalt_example_python


async def my_cor():
    return "foobar"

async def use_cpp():
    #test awaiting C++ primitives

    async for item in boost_async_example_python.test_generator():
        print("Cpp generator", item)

    print("Cpp promise gave us", await boost_async_example_python.test_promise())

    # having C++ await our python coros
    await boost_async_example_python.test_py_promise(my_cor())


asyncio.run(use_cpp())