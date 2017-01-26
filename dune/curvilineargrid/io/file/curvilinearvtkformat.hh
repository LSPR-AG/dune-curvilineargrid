#ifndef CURVILINEAR_VTK_FORMAT_HH_
#define CURVILINEAR_VTK_FORMAT_HH_



namespace Dune {

namespace CurvGrid {

/** \brief Formalism: Define XML classes, that close the bracket upon destruction. Write to XML using following statements
 *
 *  {
 *  	VTUDataArray(...)
 *  	......write data to file.....
 *  } // Destructor called here
 *
 *  */
namespace VTKFormat {

	void XMLClauseSingle(std::ofstream & file, std::string clausename, bool selfClose, bool newl = true) {
		file << "<" << clausename << (selfClose ? "/" : "") << ">";
		if (newl) { file << std::endl; }
	}

	void XMLClauseSingle(std::ofstream & file, std::string clausename, std::string clausetype, bool selfClose, bool newl = true) {
		file << "<" << clausename << " " << clausetype << (selfClose ? "/" : "") << ">";
		if (newl) { file << std::endl; }
	}

	template <class T>
	void XMLClauseSingle(std::ofstream & file, std::string clausename, std::vector<std::string> fields, std::vector<T> values, bool selfClose, bool newl = true) {
		assert(fields.size() == values.size());
		file << "<" << clausename << " ";
		for (int i = 0; i < fields.size(); i++) { file << fields[i] << "=\"" << values[i] <<"\" "; }
		file << (selfClose ? "/" : "") << ">";
		if (newl) { file << std::endl; }
	}


	struct XMLClause {
		std::ofstream & file_;
		std::string clausename_;

		XMLClause(std::ofstream & file, std::string clausename, bool newl = true) : file_(file), clausename_(clausename)
		{
			XMLClauseSingle(file, clausename, false, newl);
		}

		XMLClause(std::ofstream & file, std::string clausename, std::string clausetype, bool newl = true) : file_(file), clausename_(clausename)
		{
			XMLClauseSingle(file, clausename, clausetype, false, newl);
		}

		template <class T>
		XMLClause(std::ofstream & file, std::string clausename, std::vector<std::string> fields, std::vector<T> values, bool newl = true) : file_(file), clausename_(clausename)
		{
			XMLClauseSingle(file, clausename, fields, values, false, newl);
		}

		~XMLClause() { XMLClauseSingle(file_, "/" + clausename_, false, true); }
	};




	// [TODO] Implement binary writing, condition on format variable

	template <class BinaryType>
	class VTKArrayWriter {

	public:

		static bool isBinary(std::string format) {
			return (format == VTK_DATA_FORMAT_BINARY) || (format == VTU_DATA_FORMAT_BINARY);
		}


		template <class T>
		static void writeArray(std::ofstream & file, std::string logtext, unsigned int n, T data, std::string sep, std::string format) {
			if (isBinary(format)) {
				BinaryType binarydata = BinaryType(data);
				std::vector<char> rez(sizeof(BinaryType));
		        memcpy(rez.data(), &binarydata, sizeof(BinaryType));

				for (unsigned int i = 0; i < n; i++)  {
					if (logtext != "") { LoggingMessage::writePatience(logtext, i, n); }
					//file.write( (char*) &binarydata, sizeof(BinaryType));
					file.write (rez.data(), rez.size());
				}
			} else {
				for (unsigned int i = 0; i < n; i++)  {
					if (logtext != "") { LoggingMessage::writePatience(logtext, i, n); }
					file << data << sep;
				}
			}
		}


		template <class T>
		static void writeArray(std::ofstream & file, std::string logtext, const std::vector<T> & dataVec, std::string sep, std::string format)
		{
			if (isBinary(format)) {
				for (unsigned int i = 0; i < dataVec.size(); i++)  {
					if (logtext != "") { LoggingMessage::writePatience(logtext, i, dataVec.size()); }
					BinaryType binarydata = BinaryType(dataVec[i]);
					std::vector<char> rez(sizeof(BinaryType));
			        memcpy(rez.data(), &binarydata, sizeof(BinaryType));
			        file.write (rez.data(), rez.size());

					//file.write( (char*) &binarydata, sizeof(BinaryType));
				}
			} else {
				for (unsigned int i = 0; i < dataVec.size(); i++)  {
					if (logtext != "") { LoggingMessage::writePatience(logtext, i, dataVec.size()); }
					file << dataVec[i] << sep;
				}
			}

		}


		template <int dim, class ctype>
		static void writeCoordinateArray(std::ofstream & file, std::string logtext, const std::vector<Dune::FieldVector<ctype,dim>> & dataVec, std::string sep, std::string format)
		{
			if (isBinary(format)) {
				for (unsigned int i = 0; i < dataVec.size(); i++)  {
					if (logtext != "") { LoggingMessage::writePatience(logtext, i, dataVec.size()); }
					for (unsigned int d = 0; d < dim; d++) {
						BinaryType binarydata = BinaryType(dataVec[i][d]);
						std::vector<char> rez(sizeof(BinaryType));
				        memcpy(rez.data(), &binarydata, sizeof(BinaryType));
				        file.write (rez.data(), rez.size());
						//file.write( (char*) &binarydata, sizeof(BinaryType));
					}
				}
			} else {
				for (unsigned int i = 0; i < dataVec.size(); i++)  {
					if (logtext != "") { LoggingMessage::writePatience(logtext, i, dataVec.size()); }
					for (unsigned int d = 0; d < dim; d++) {
						file << dataVec[i][d] << " ";
					}
					file << sep;
				}
			}
		}
	};



} // Namespace VTKFormat

} // Namespace CurvGrid

} // Namespace Dune

#endif // CURVILINEAR_VTK_FORMAT_HH_
