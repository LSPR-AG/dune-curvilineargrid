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

	void XMLClauseSingle(std::ofstream & file, std::string clausename, bool selfClose) {
		file << "<" << clausename << (selfClose ? "/" : "") << ">" << std::endl;
	}

	void XMLClauseSingle(std::ofstream & file, std::string clausename, std::string clausetype, bool selfClose) {
		file << "<" << clausename << " " << clausetype << (selfClose ? "/" : "") << ">" << std::endl;
	}

	template <class T>
	void XMLClauseSingle(std::ofstream & file, std::string clausename, std::vector<std::string> fields, std::vector<T> values, bool selfClose) {
		assert(fields.size() == values.size());
		file << "<" << clausename << " ";
		for (int i = 0; i < fields.size(); i++) { file << fields[i] << "=\"" << values[i] <<"\" "; }
		file << (selfClose ? "/" : "") << ">" << std::endl;
	}


	struct XMLClause {
		std::ofstream & file_;
		std::string clausename_;

		XMLClause(std::ofstream & file, std::string clausename) : file_(file), clausename_(clausename)
		{
			XMLClauseSingle(file, clausename, false);
		}

		XMLClause(std::ofstream & file, std::string clausename, std::string clausetype) : file_(file), clausename_(clausename)
		{
			XMLClauseSingle(file, clausename, clausetype, false);
		}

		template <class T>
		XMLClause(std::ofstream & file, std::string clausename, std::vector<std::string> fields, std::vector<T> values) : file_(file), clausename_(clausename)
		{
			XMLClauseSingle(file, clausename, fields, values, false);
		}

		~XMLClause() { XMLClauseSingle(file_, "/" + clausename_, false); }
	};


	// [TODO] Implement binary writing, condition on format variable
	template <class T>
	void writeArray(std::ofstream & file, std::string logtext, unsigned int n, T data, std::string sep, std::string format) {
		for (unsigned int i = 0; i < n; i++)  {
			if (logtext != "") {LoggingMessage::writePatience(logtext, i, n);}
			file << data << sep;
		}
	}

	template <class T>
	void writeArray(std::ofstream & file, std::string logtext, const std::vector<T> & dataVec, std::string sep, std::string format)
	{
		for (unsigned int i = 0; i < dataVec.size(); i++)  {
			if (logtext != "") { LoggingMessage::writePatience(logtext, i, dataVec.size()); }
			file << dataVec[i] << sep;
		}
	}

	template <int dim, class ctype>
	void writeCoordinateArray(std::ofstream & file, std::string logtext, const std::vector<Dune::FieldVector<ctype,dim>> & dataVec, std::string sep, std::string format)
	{
		for (unsigned int i = 0; i < dataVec.size(); i++)  {
			if (logtext != "") { LoggingMessage::writePatience(logtext, i, dataVec.size()); }
			for (unsigned int d = 0; d < dim; d++) {
				file << dataVec[i][d] << " ";
			}
			file << sep;
		}
	}

} // Namespace VTKFormat

} // Namespace CurvGrid

} // Namespace Dune

#endif // CURVILINEAR_VTK_FORMAT_HH_
