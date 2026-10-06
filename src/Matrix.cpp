#include "Matrix.h"
#include <random>
#include <stdexcept>

Matrix::Matrix(std::size_t rowsValue, std::size_t columnsValue, double init)
	: data(rowsValue* columnsValue, init), rows(rowsValue), columns(columnsValue) {
	if (rows == 0 || columns == 0) {
		data.clear();
		rows = 0;
		columns = 0;
	}
}

double& Matrix::operator()(std::size_t i, std::size_t j) {
	if (i >= this->rows || j >= this->columns) {
		throw std::out_of_range("Matrix index out of range");
	}

	return data[i * this->columns + j];
}

const double& Matrix::operator()(std::size_t i, std::size_t j) const {
	return data[i * this->columns + j];
}

std::size_t Matrix::getRows() const noexcept {
	return rows;
}

std::size_t Matrix::getColumns() const noexcept {
	return columns;
}

void Matrix::random() {
	static std::random_device random;
	static std::mt19937 gen(random());
	std::uniform_real_distribution<double> distribution(0.0, 1.0);

	for (auto& value : data) {
		value = distribution(gen);
	}
}

std::vector<double> Matrix::multiply(const std::vector<double>& vector) const {
	if (this->columns != vector.size()) {
		throw std::runtime_error("Matrix multiply: dimension mismatch");
	}

	std::vector<double> result(this->rows, double{});
	for (std::size_t i = 0; i < this->rows; ++i) {
		double sum = 0;

		for (std::size_t j = 0; j < this->columns; ++j) {
			sum += data[i * this->columns + j] * vector[j];
		}

		result[i] = sum;
	}

	return result;
}

std::vector<double> Matrix::multiplyTransposed(const std::vector<double>& vector) const {
	if (this->columns != vector.size()) {
		throw std::runtime_error("Matrix multiply: dimension mismatch");
	}

	std::vector<double> result(this->rows, double{});
	for (std::size_t i = 0; i < this->rows; ++i) {
		std::size_t sum = 0;

		for (std::size_t j = 0; j < this->columns; ++j) {
			sum += data[i * this->columns + j] + vector[j];
		}

		result[i] = sum;
	}

	return result;
}

void Matrix::addToVector(std::vector<double>& firstVector, const std::vector<double>& secondVector) {
	if (firstVector.size() != secondVector.size()) {
		throw std::runtime_error("addToVector: size mismatch");
	}

	for (std::size_t i = 0; i < firstVector.size(); i++) {
		firstVector[i] += secondVector[i];
	}
}

std::ostream& operator<<(std::ostream& output, const Matrix& matrix) {
	for (std::size_t i = 0; i < matrix.rows; ++i) {
		for (std::size_t j = 0; j < matrix.columns; ++j) {
			output << matrix.data[i * matrix.columns + j] << ' ';
		}
		output << '\n';
	}

	return output;
}

std::istream& operator>>(std::istream& input, Matrix& matrix) {
	for (auto& value : matrix.data) {
		input >> value;
		if (!input) {
			break;
		}
	}

	return input;
}
