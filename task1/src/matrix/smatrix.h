#ifndef SMATRIX_H
#define SMATRIX_H

#define USE_DOUBLE 1
#define USE_COMPLEX 1

#include "mpi.h"
#include <ostream>
#include <complex>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

#ifdef USE_DOUBLE
	typedef double SReal;
	#define pRgesvd pdgesvd_
	#define pCgesvd pzgesvd_
#else
	typedef float SReal;
	#define pRgesvd psgesvd_
	#define pCgesvd pcgesvd_
#endif

#ifdef USE_COMPLEX
	#define pNgesvd pCgesvd
	typedef complex<SReal> SType;
#else
	#define pNgesvd pRgesvd
	typedef SReal SType;
#endif

class SMatrix {
public:
	/** \brief The constructor for SMatrix.
	* \param nRows Number of rows in the matrix.
	* \param nCols Number of columns in the matrix.
	* \param nprocRows The numer of rows in the processes grid. \c 0 means auto guess.
	* \param nprocCols The number of columns in the processes grid. \c 0 means auto guess.
	* 
	* The total number of used processes is \a nprocRows * \a nprocCols.
	* The number of used processes should be less than or equal to the number of total processes that the program has.
	* \sa init(), ~SMatrix()
	*/
	SMatrix(int nRows, int nCols, int nProcRows=0, int nProcCols=0);
	
	/** \brief The copy constructor for SMatrix.
	* \param matrix The matrix to copy.
	* 
	* Copies the given SMatrix \p matrix into a new one.
	* The new matrix is created in the same context as the original one.
	* \sa init(), ~SMatrix()
	*/
	SMatrix(const SMatrix& matrix);

	/** \brief The same context constructor for SMatrix.
	* \param matrix The parent matrix.
	* \param nRows Number of rows in the new matrix.
	* \param nCols Number of columns in the new matrix.
	* 
	* Creates a new matrix in the same context as the parent one.
	* The new matrix will have the same \a nprocRows and \a nprocCols as the parent one.
	* \sa init(), ~SMatrix()
	*/
	SMatrix(const SMatrix& matrix, int nRows, int nCols);


	/** \brief The destructor for SMatrix.
	*
	* This method destucts the SMatrix object.
	* \sa exit(), SMatrix()
	*/
	~SMatrix();

	/** \brief Initalize the processing context.
	* \param[out] rank A pointer to an integer where the rank of the current process is written.
	* \param[out] nprocs A pointer to an integer where the number of all processes is written.
	* \param blockSize Blocking size. Rows and columns will be distributed among processes in portions of that size.
	* 
	* This method initializes the processing context.
	* This method should be called before any other uses of SMatrix.
	* \return \c True if the current process is root, \c False otherwise.
	* \sa exit(), SMatrix()
	*/
	static bool init(int *rank=NULL, int *nprocs=NULL, int blockSize=64);

	/** \brief Destroy the processing context.
	*
	* This method destroys the processing context.
	* This method should be called after any other uses of SMatrix.
	* Using SMatrix after calling this method is not permitted.
	* \sa init(), ~SMatrix()
	*/
	static void exit();

	/** \brief Wait for all processes.
	 *
	 * Produces a barrier that waits for all processes.
	 */
	static void barrier();

	/** \brief Dump internal information.
	 *
	 * A debug function that makes all processes dump their internal information.
	 */
	void dumpInfo() const;

	/** \brief Dump all data blocks one by one.
	 *
	 * A debug function that makes all processes dump their parts of the table one by one.
	 */
	void dump() const;

	/** \brief Formatted output.
	 * \param stream The stream in which the matrix should be written.
	 * 
	 * Output the formatted matrix to the \p stream.
	 */
	void operator>>(ostream& stream) const;
	

	/** \brief Fill the specified block with given data.
	 * \param procRow The number of the row in the processes grid.
	 * \param procCol The number of the column in the processes grid.
	 * \param[in] data The pointer to the data array.
	 * \param size How many items should be copied from \p data. \c -1 means copying everything to the end of the block.
	 * \param offset The offset in the block from which the elements should be copied to the block.
	 * 
	 * Fills the block at [\p procRow][\p procCol] with \p size items from \p data, starting from \p offset.
	 * \warning Does not produce a barrier.
	 * \sa set()
	 */
	void setBlock(int procRow, int procCol, const SType *data, int size = -1, int offset = 0);

	/** \brief Copy the data to the matrix.
	 * \param[in] data The pointer to the data array.
	 * \param rowStart The number of the first row of the submatrix.
	 * \param colStart The number of the first column of the submatrix.
	 * \param rowCount The number of the rows to copy. \c 0 means "to the end of the matrix".
	 * \param colCount The number of the columns to copy. \c 0 means "to the end of the matrix".
	 * 
	 * Fills the submatrix at [\p rowStart][\p colStart] with size [\p rowCount][\p colCount] with values from \p data.
	 * \warning Does not produce a barrier.
	 * \sa setBlock()
	 */
	void set(const SType *data, int rowStart = 0, int colStart = 0, int rowCount = 0, int colCount = 0);


	/** \brief Fill the matrix with a given number.
	 * \param value The value with which the matrix will be filled.
	 * 
	 * Fills the whole matrix with \p value.
	 * \sa setZero(), populate()
	 * \sa operator=(SType value)
	 */
	void fill(SType value);

	/** \brief Populate the matrix using a given function.
	 * \param[in] func The function that returns the values by their indexes.
	 * 
	 * Populates the whole matrix using the function \p func.
	 * \sa fill()
	 */
	void populate(SType (*func)(int, int));

	/** \brief Zero out the matrix.
	 *
	 * Fills the whole matrix with zeros.
	 * \sa setIdentity(), fill()
	 */
	void setZero();

	/** \brief Set matrix to identity.
	 *
	 * Sets the matrix to an identity matrix of the same size, deleting previous data.
	 * \sa setZero(), fill(), populate()
	 */
	void setIdentity();


	/** \brief Copy operator.
	 * 
	 * Copy given matrix \p matrix into current one.
	 * The context is shared after the copy operation.
	 * \return The reference to the matrix.
	 * \sa SMatrix(const SMatrix&)
	 * \sa operator==(const SMatrix&)
	 */
	SMatrix& operator =(const SMatrix& matrix);

	/** \brief Compare operator.
	 * \param matrix The matrix to compare with.
	 * 
	 * Compare the matrix with a given matrix \p matrix.
	 * \return \c True, if the matrix is equal to \p matrix, \c False otherwise.
	 * \warning Compares memory. You probably want fuzzyCompare() on floating point matrices.
	 * \sa operator=(const SMatrix&)
	 */
	bool operator ==(const SMatrix& matrix);

	/** \brief Fill the matrix with a given number.
	 * \param value The value with which the matrix will be filled.
	 * 
	 * This is an alias for fill(SType value).
	 * \return The reference to the matrix.
	 * \sa operator*=(SType value), operator/=(SType value), operator+=(SType value), operator-=(SType value)
	 * \sa fill()
	 */
	SMatrix& operator =(SType value);

	/** \brief Multiply the matrix by a number.
	 * \param value The multiplier.
	 * 
	 * Multiplies the whole matrix by the given number \p value.
	 * \return The reference to the matrix.
	 * \sa operator/=(SType value), operator+=(SType value). operator-=(SType value), operator=(SType value)
	 */
	SMatrix& operator *=(SType value);

	/** \brief Divide the matrix by a number.
	 * \param value The divisor.
	 * 
	 * Divides the whole matrix by the given number \p value.
	 * \return The reference to the matrix.
	 * \sa operator*=(SType value), operator+=(SType value), operator-=(SType value), operator=(SType value)
	 */
	SMatrix& operator /=(SType value);

	/** \brief Add a number to all elements.
	 * \param value The number to add to all elements.
	 * 
	 * Adds the given number \p value to all elements of the matrix.
	 * \return The reference to the matrix.
	 * \sa operator-=(SType value), operator*=(SType value), operator/=(SType value), operator=(SType value)
	 */
	SMatrix& operator +=(SType value);

	/** \brief Subtract a number from all elements.
	 * \param value The number to subtract from all elements.
	 * 
	 * Subtracts the given number \p value from all elements of the matrix.
	 * \return The reference to the matrix.
	 * \sa operator+=(SType value), operator*=(SType value), operator/=(SType value), operator=(SType value)
	 */
	SMatrix& operator -=(SType value);

	/** \brief Calculate singular value decomposition.
	 * \param[out] SSize The size of the output array of singular numbers.
	 * \param[out] U The pointer to the left unitary matrix.
	 * \param[out] VT The pointer to the right unitary matrix, transposed.
	 * 
	 * Calculates singular value decomposition and returns all parts of it.
	 * \return The array of singular values.
	 */
	SReal* calculateEIGENVAL(int *SSize, SMatrix **Z = NULL);
	void transpconj(bool isRoot);
	void mul(char *ch, char *ch_B, SMatrix *B, SMatrix **C);
	void DiagToMat(int *SSize, SType *D);
	void readf(MPI_File thefile);
	void wrbuf(SType* buf);
	void setij(int i, int j);

private:
	void create(int icon = -2);
	void destroy();
	
	static SType fillerIdentity(int i, int j);
	int* getDesc();
	int* getDesc(int rows, int cols, int procRows);

	static int rank, nProcs, blockSize, root;
	int icon;
	int nRows, nCols, nProcRows, nProcCols,
		myProcRow, myProcCol,
		myProcRowsOffset, myProcColsOffset,
		myProcRows, myProcCols, myProcSize;
	bool myProcOk;
	SType *data;
};

ostream& operator<<(ostream& out, const SMatrix& matrix);

#endif
