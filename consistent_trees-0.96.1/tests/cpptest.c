#undef A
#define A 1
#undef A
#define A foo

int A(int i) {
  return i;
}

int main(void) {
  return A(1);
}
