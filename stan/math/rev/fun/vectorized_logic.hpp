#include <vector>

using namespace std;

// Vector-vector equality
template <typename T>
vector<int> logical_eq(const vector<T>& vector1, const vector<T>& vector2) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] == vector2[i]);
      }
    return result;
}

// Vector-vector less than
template <typename T>
vector<int> logical_lt(const vector<T>& vector1, const vector<T>& vector2) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] < vector2[i]);
      }
    return result;
}

// Vector-vector greater than
template <typename T>
vector<int> logical_gt(const vector<T>& vector1, const vector<T>& vector2) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] > vector2[i]);
      }
    return result;
}

// Vector-vector less than or equal to
template <typename T>
vector<int> logical_leq(const vector<T>& vector1, const vector<T>& vector2) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] <= vector2[i]);
      }
    return result;
}

// Vector-vector greater than or equal to
template <typename T>
vector<int> logical_geq(const vector<T>& vector1, const vector<T>& vector2) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] >= vector2[i]);
      }
    return result;
}
// Vector-scalar equality
template <typename T>
vector<int> logical_eq(const vector<T>& vector1, T scalar) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] == scalar);  
      }
    return result;
}

// Vector-scalar less than
template <typename T>
vector<int> logical_lt(const vector<T>& vector1, T scalar) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] < scalar);  
      }
    return result;
}

// Vector-scalar greater than
template <typename T>
vector<int> logical_gt(const vector<T>& vector1, T scalar) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] > scalar);  
      }
    return result;
}

// Vector-scalar less than or equal
template <typename T>
vector<int> logical_leq(const vector<T>& vector1, T scalar) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] <= scalar);  
      }
    return result;
}

// Vector-scalar greater than or equal
template <typename T>
vector<int> logical_geq(const vector<T>& vector1, T scalar) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] >= scalar);  
      }
    return result;
}

// Vector element-wise negation
template <typename T>
vector<int> logical_not(const vector<T>& vector1) {
    vector<int> result(vector1.size());
    for (size_t i = 0; i < vector1.size(); i++) 
      {
        result[i] = (vector1[i] == 0);
      }
    return result;
}

// Vector size-wise OR
template <typename T>
int any(const vector<T>& vector1) {
    for (size_t i = 0; i < vector1.size(); i++) {
        if (vector1[i] != 0) 
        {
          return 1;
        };
    }
    return 0;
}

// Vector size-wise AND
template <typename T>
int all(const vector<T>& vector1) {
    for (size_t i = 0; i < vector1.size(); i++) {
        if (vector1[i] == 0) 
        {
          return 0;
        }
    }
    return 1;
}
