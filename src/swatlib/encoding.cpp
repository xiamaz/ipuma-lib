#include "encoding.h"
#include <stdexcept>

namespace swatlib {
	DataType strToDataType(std::string s) {
		if (s == "na") {
			return DataType::nucleicAcid;
		} else if (s == "aa") {
			return DataType::aminoAcid;
		} else if (s == "str") {
			return DataType::string;
		} else {
			throw std::runtime_error("Invalid dataType: " + s);
		}
	}

	std::string dataTypeToStr(DataType d) {
		switch (d) {
			case DataType::aminoAcid:
				return "aa";
				break;
			case DataType::nucleicAcid:
				return "na";
				break;
			case DataType::string:
				return "str";
				break;
		}
		throw std::runtime_error("Invalid DataType");
	}

               Encoding::Encoding(std::map<char, uint8_t> mapping) : code(mapping) {
		for (auto& [key, val] : code) {
			inverse[val] = key;
		}
	}
	uint8_t Encoding::encode(char c) {
		return code[c];
	}
	char Encoding::decode(uint8_t i) {
		return inverse[i];
	}

	std::vector<uint8_t> Encoding::encode(const std::string& s) {
		std::vector<uint8_t> v(s.size());
		for (int i = 0; i < s.size(); ++i) {
			v[i] = this->encode(s[i]);
		}
		return v;
	}

	std::vector<std::vector<uint8_t>> Encoding::encode(const std::vector<std::string>& ss) {
		std::vector<std::vector<uint8_t>> vv(ss.size());
		for (int i = 0; i < ss.size(); ++i) {
			vv[i] = encode(ss[i]);
		}
		return vv;
	}


	std::string Encoding::decode(const std::vector<uint8_t>& v) {
		std::string s(v.size(), ' ');
		for (int i = 0; i < v.size(); ++i) {
			s[i] = this->decode(v[i]);
		}
		return s;
	}

	std::vector<std::string> Encoding::decode(const std::vector<std::vector<uint8_t>>& vv) {
		std::vector<std::string> ss(vv.size());
		for (int i = 0; i < vv.size(); ++i) {
			ss[i] = decode(vv[i]);
		}
		return ss;
	}

	uint8_t Encoding::get_terminator() {
		return code['\0'];
	}

	NAEncoding::NAEncoding() : Encoding::Encoding(NA_CODE) {}
	AAEncoding::AAEncoding() : Encoding::Encoding(AA_CODE) {}

	std::map<char, uint8_t> createCharCode() {
		std::map<char, uint8_t> m;
		for (int i = 0; i < 256; ++i) {
			m[static_cast<char>(i)] = i;
		}
		return m;
	}

	CharEncoding::CharEncoding() : Encoding::Encoding(createCharCode()) {}

	Encoding getEncoder(DataType t) {
		switch(t) {
			case DataType::aminoAcid:
			return AAEncoding();
			break;
			case DataType::nucleicAcid:
			return NAEncoding();
			break;
			case DataType::string:
			return CharEncoding();
			break;
		}
		throw std::runtime_error("Invalid DataType");
	}

}