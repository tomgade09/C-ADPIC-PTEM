#include "utils/serializationHelpers.h"

namespace utils
{
	namespace fileIO
	{
		namespace serialize
		{
			size_t readSizetLength(ifstream& in)
			{
				vector<char> size(sizeof(size_t));
				size_t ret{ 0 };

				in.read(size.data(), sizeof(size_t));
				ret = *(reinterpret_cast<size_t*>(size.data()));

				return ret;
			}

			// ================ serialize functions ================ //
			stringbuf serializeString(const string& str)
			{
				stringbuf sb;
				ostream out(&sb);

				size_t size{ str.size() };
				out.write(reinterpret_cast<char*>(&size), sizeof(size_t));

				out << str;

				return sb;
			}

			stringbuf serializeDoubleVector(const vector<double>& vec)
			{
				stringbuf sb;
				ostream out(&sb);

				size_t size{ vec.size() };
				out.write(reinterpret_cast<char*>(&size), sizeof(size_t));

				for (auto& elem : vec)
					out.write(reinterpret_cast<const char*>(&elem), sizeof(double));

				return sb;
			}

			stringbuf serializeStringVector(const vector<string>& vec)
			{
				stringbuf sb;
				ostream out(&sb);

				size_t size{ vec.size() };
				out.write(reinterpret_cast<char*>(&size), sizeof(size_t));

				for (auto& str : vec)
					out << serializeString(str).str();

				return sb;
			}

			// ================ deserialize functions ================ //
			string deserializeString(ifstream& istr)
			{
				size_t strlen{ readSizetLength(istr) };

				vector<char> strchar(strlen);
				istr.read(strchar.data(), strlen);

				string ret(strchar.data(), strlen);

				return ret;
			}

			vector<double> deserializeDoubleVector(ifstream& istr)
			{
				size_t veclen{ readSizetLength(istr) };
				vector<double> ret;

				for (size_t elem = 0; elem < veclen; elem++)
				{
					vector<char> dblchar(sizeof(double));
					istr.read(dblchar.data(), sizeof(double));
					ret.push_back(*(reinterpret_cast<double*>(dblchar.data())));
				}

				return ret;
			}

			vector<string> deserializeStringVector(ifstream& istr)
			{
				size_t vecsize{ readSizetLength(istr) };
				vector<string> ret;

				for (size_t elem = 0; elem < vecsize; elem++)
					ret.push_back(deserializeString(istr));

				return ret;
			}
		}
	}
}