
namespace GP {

	float min(float x, float y) {
		return (x < y) ? x : y;
	}

	void convmtx(arma::vec& kern, int nsamp, arma::mat& convmat) {
		// convolution matrix
		int k, ind;
		int ncol = kern.size();
		int nrow = ncol + nsamp - 1;
		convmat.zeros(nrow, nrow);
		k = 0;

		for (int i = 0; i < nrow; i++) {
			ind = ncol + k;

			if (ind > nrow)
			{
				ind = nrow;
			}

			for (int ii = 0 + k; ii < ind; ii++)
			{
				convmat(ii, i) = kern[ii - k];
			}

			k = k + 1;
		}
	}

	void convmtx2(arma::vec& kern, int nsamp, arma::mat& convmat) {
		// convolution matrix

		arma::mat temp;
		convmat.zeros(nsamp, nsamp);

		// calls convolution matrix
		convmtx(kern, nsamp, temp);
		convmat = temp(arma::span(0, nsamp - 1), arma::span(0, nsamp - 1));

	}

	void dervmtx(arma::mat& dermat, int nsamp) {
		// derivative matrix

		dermat.zeros(nsamp, nsamp);

		for (int i = 0; i < nsamp - 2; i++) {
			dermat(i, i) = -1;
			dermat(i, i + 1) = 1;
		}
		dermat(nsamp - 1, nsamp - 1) = -1;


	}

	void calc_AI(arma::vec& vec1, arma::vec& vec2, arma::vec& AI) {

		int ns1 = vec1.size();
		int ns2 = vec2.size();
		AI.zeros(ns1);

		if (ns1 != ns2) {
			std::cout << "Unequal vector length" << std::endl;
			//		exit;
		}

		for (int i = 0; i < ns1; i++) {
			AI[i] = vec1(i)*vec2(i);
		}

	}

	void calc_RC(arma::vec& AI, arma::vec& RCF) {

		int ns = AI.size();
		RCF.zeros(ns);



		for (int i = 0; i < ns-1; i++) {
			RCF[i] = (AI(i + 1) - AI(i)) / (AI(i + 1) + AI(i));
		}

		RCF(ns - 1) = RCF(ns - 2);

	}

	void conv(arma::vec& signal, arma::vec& kernel, arma::vec& filtered) {
		int j, k;
		int datasize = signal.size();
		int kernelsize = kernel.size(); 


		// last part
		for (int i = kernelsize; i < datasize; i++) {
			filtered[i] = 0;
			j = i;

			for (int k = 0; i < kernelsize; k++) {
				filtered[i] = signal[i] * signal[j] * kernel[k];
				j--;
			}
		}

		// first part

		for (int i = 0; i < kernelsize; i++) {
			filtered[i] = 0;
			j = i;
			k = 1;
			while (j > 0) {
				filtered[i] = signal[i] * signal[j] * kernel[k];
				j--;
				k++;
			}
		}
	}








	//---------------Stats-------------------------


	float corr(arma::vec vec1, arma::vec vec2) {

		float r, corr, meanvec1, meanvec2;
		float tmp1, tmp2, tmp12;
		int lag;

		tmp1 = arma::dot(vec1 - mean(vec1), vec1 - mean(vec1));
		tmp2 = arma::dot(vec2 - mean(vec2), vec2 - mean(vec2));
		tmp12 = arma::dot(vec1 - mean(vec1), vec2 - mean(vec2));

		r = tmp12 / (std::sqrt(tmp1)*std::sqrt(tmp2));

		return r;

	}



	float sum(arma::vec& vec1, int nsamp) {


		float summ = 0;
		for (int i = 0; i < nsamp; i++)
		{
			summ = summ + vec1[i];
		}

		return summ;
	}


	float sum(arma::vec& vec1) {

		int nsamp = vec1.size();

		float summ = 0;
		for (int i = 0; i < nsamp; i++)
		{
			summ = summ + vec1[i];
		}

		return summ;
	}


	float avg(arma::vec& vec1) {

		float tmp;
		int nsamp = vec1.size();
		tmp = sum(vec1);
		return (tmp / nsamp);

	}


	float avg(arma::vec& vec1, int nsamp) {
		float tmp;
		tmp = sum(vec1, nsamp);
		return (tmp / nsamp);
	}


	void span(arma::mat& mtx, int m_col, int m_begin, int m_end, arma::vec& vec1, int v_begin, int v_end)
	{
		int m_nsamp = m_end - m_begin;
		int v_nsamp = v_end - v_begin;
		if (m_nsamp != v_nsamp)
		{
			std::cout << "size are different , please check" << std::endl;
			std::cout << "matrix row size  = " << m_nsamp << std::endl;
			std::cout << "vector row size  = " << v_nsamp << std::endl;
			exit(0);
		}

		for (int i = 0; i < m_nsamp; i++)
		{

			mtx(m_col, m_begin + i) = vec1(v_begin + i);
		}
	}

	void indM2V(arma::mat& mtx, int m_col, int m_begin, int m_end,
		arma::vec& vec1, arma::ivec indx, int v_begin, int v_end)
	{
		int m_nsamp = m_end - m_begin;
		int v_nsamp = v_end - v_begin;
		float tmpval, tmp_ind;


		for (int i = 0; i < m_nsamp; i++)
		{
			tmp_ind = indx[v_begin + i];
			tmpval = vec1[tmp_ind];
			mtx(m_begin + i, m_col) = tmpval;

			//		std::cout << "value of back-transformed AI = " << tmpval << "at index = " << tmp_ind << std::endl;
		}


	}

	void indV2M(arma::mat& mtx, int m_col, int m_begin, int m_end,
		arma::vec& vec1, arma::ivec indx, int v_begin, int v_end)
	{
		int m_nsamp = m_end - m_begin ;
		int v_nsamp = v_end - v_begin ;
		float tmpval, tmp_ind;


		for (int i = 0; i < m_nsamp; i++)
		{
			tmp_ind = indx[v_begin + i];
			tmpval = vec1[tmp_ind];
			mtx(m_begin + i, m_col) = tmpval;

	//		std::cout << "value of back-transformed AI = " << tmpval << "at index = " << tmp_ind << std::endl;
		}


	}



	void indV2V(arma::vec& vec_out, int m_begin, int m_end,
		arma::vec& vec1, arma::ivec indx, int v_begin, int v_end)
	{
		int m_nsamp = m_end - m_begin;
		int v_nsamp = v_end - v_begin;
		float tmpval, tmp_ind;

		if (m_nsamp != v_nsamp)
		{
			std::cout << "size are different , please check" << std::endl;
			std::cout << "vec_out  size  = " << m_nsamp << std::endl;
			std::cout << "vec_ind  size  = " << v_nsamp << std::endl;
			exit(0);
		}

		for (int i = 0; i < m_nsamp; i++)
		{
			tmp_ind = indx[v_begin + i];
			tmpval = vec1[tmp_ind];
			vec_out(m_begin + i) = tmpval;
		}


	}




	float*  arma2array(arma::vec& vec1) {

		float*arr = NULL; // note NULL = #define = o 
		int ns = vec1.size();

		arr = new float[ns];

		for (int i = 0; i < ns; i++)
		{
			arr[i] = vec1[i];
		}

		return arr;

	}
}