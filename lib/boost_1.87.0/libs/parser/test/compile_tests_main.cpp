// Copyright (C) 2018 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
void compile_attribute();
void compile_seq_attribute();
void compile_or_attribute();
void compile_combining_groups();
void compile_all_t();

int main()
{
    compile_attribute();
    compile_seq_attribute();
    compile_or_attribute();
    compile_combining_groups();
    compile_all_t();
}
