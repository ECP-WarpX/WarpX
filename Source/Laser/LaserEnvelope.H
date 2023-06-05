#include <complex>

int main() {
    std::complex<double> A_envelope_potential_vector;

    // Accessing and modifying the real and imaginary parts
    A_envelope_potential_vector.real(0.0);
    A_envelope_potential_vector.imag(0.0);

    // Printing the real and imaginary parts
    std::cout << "Real part: " << myComplexField.real() << std::endl;
    std::cout << "Imaginary part: " << myComplexField.imag() << std::endl;

    return 0;
}
