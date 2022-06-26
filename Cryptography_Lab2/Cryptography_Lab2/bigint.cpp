#include "bigint.h"

bigint::bigint()
{
	this->_is_negative = false;
}

bigint::bigint(std::string str)
{
	if (str.length() == 0)
	{
		this->_is_negative = false;
	}
	else
	{
		if (str[0] == '-')
		{
			str = str.substr(1);
			this->_is_negative = true;
		}
		else
		{
			this->_is_negative = false;
		}
		copyStr = str;
		for (long long i = str.length(); i > 0; i -= 9)
		{
			if (i < 9)
				this->_digits.push_back(atoi(str.substr(0, i).c_str()));
			else
				this->_digits.push_back(atoi(str.substr(i - 9, 9).c_str()));
		}

		this->_remove_leading_zeros();
	}
}

string bigint::toStr() // перевод в string
{
	vector<int> destination;
	stringstream result;
	for (auto it = _digits.rbegin(); it != _digits.rend(); ++it)
	{
		destination.push_back(*it);
	}
	copy(destination.begin(), destination.end(), ostream_iterator<int>(result, ""));

	return result.str().c_str();;
}

void bigint::_remove_leading_zeros() // ігноруємо 0
{
	while (this->_digits.size() > 1 && this->_digits.back() == 0)
	{
		this->_digits.pop_back();
	}

	if (this->_digits.size() == 1 && this->_digits[0] == 0) this->_is_negative = false;
}

std::ostream& operator <<(std::ostream& os, const bigint& bi) // друкує число в потік виведення
{
	if (bi._digits.empty()) os << 0;
	else {
		if (bi._is_negative) os << '-';
		os << bi._digits.back();
		char old_fill = os.fill('0');
		for (long long i = static_cast<long long>(bi._digits.size()) - 2; i >= 0; --i)
			os << std::setw(9) << bi._digits[i];
		os.fill(old_fill);
	}

	return os;
}

bool operator ==(const bigint& left, const bigint& right) // оператор порівняння
{
	if (left._is_negative != right._is_negative) return false;
	if (left._digits.empty())
	{
		if (right._digits.empty() || (right._digits.size() == 1 && right._digits[0] == 0)) return true;
		else return false;
	}

	if (right._digits.empty())
	{
		if (left._digits.size() == 1 && left._digits[0] == 0) return true;
		else return false;
	}

	if (left._digits.size() != right._digits.size()) return false;
	for (size_t i = 0; i < left._digits.size(); ++i) if (left._digits[i] != right._digits[i]) return false;

	return true;
}

const bigint bigint::operator +() const // унарний +
{
	return bigint(*this);
}

const bigint bigint::operator -() const // унарний -
{
	bigint copy(*this);
	copy._is_negative = !copy._is_negative;
	return copy;
}

bool operator <(const bigint& left, const bigint& right) // оператор менше
{
	if (left == right) return false;
	if (left._is_negative)
	{
		if (right._is_negative) return ((-right) < (-left));
		else return true;
	}
	else if (right._is_negative) return false;
	else
	{
		if (left._digits.size() != right._digits.size())
		{
			return left._digits.size() < right._digits.size();
		}
		else
		{
			for (long long i = left._digits.size() - 1; i >= 0; --i)
			{
				if (left._digits[i] != right._digits[i]) return left._digits[i] < right._digits[i];
			}
			return false;
		}
	}
}

bool operator !=(const bigint& left, const bigint& right)// порівнює два числа на нерівність
{
	return !(left == right);
}

bool operator <=(const bigint& left, const bigint& right) // переіврка <=
{
	return (left < right || left == right);
}

bool operator >(const bigint& left, const bigint& right) // перевірка >
{
	return !(left <= right);
}

bool operator >=(const bigint& left, const bigint& right) // перевірка >=
{
	return !(left < right);
}

const bigint operator +(bigint left, const bigint& right) // додавання двух чисел
{
	if (left._is_negative)
	{
		if (right._is_negative) return -(-left + (-right));
		else return right - (-left);
	}
	else if (right._is_negative) return left - (-right);
	int carry = 0;
	size_t max = 0;
	if (left._digits.size() >= right._digits.size())
	{
		max = left._digits.size();
	}
	else
	{
		max = right._digits.size();
	}
	for (size_t i = 0; i < max || carry != 0; ++i)
	{
		if (i == left._digits.size()) left._digits.push_back(0);
		left._digits[i] += carry + (i < right._digits.size() ? right._digits[i] : 0);
		carry = left._digits[i] >= bigint::BASE;
		if (carry != 0) left._digits[i] -= bigint::BASE;
	}

	return left;
}

bigint& bigint::operator +=(const bigint& value)
{
	return *this = (*this + value);
}

const bigint bigint::operator++() // префіксний інкремент
{
	return (*this += 1);
}

bigint::operator std::string() const // перетворює число в string
{
	std::stringstream ss;
	ss << *this;
	return ss.str();
}

bigint::bigint(signed char c) {
	if (c < 0) this->_is_negative = true;
	else this->_is_negative = false;
	this->_digits.push_back(std::abs(c));
}

bigint::bigint(unsigned char c)
{
	this->_is_negative = false;
	this->_digits.push_back(c);
}

bigint::bigint(signed short s)
{
	if (s < 0) this->_is_negative = true;
	else this->_is_negative = false;
	this->_digits.push_back(std::abs(s));
}

bigint::bigint(unsigned short s)
{
	this->_is_negative = false;
	this->_digits.push_back(s);
}

bigint::bigint(signed int i)
{
	if (i < 0) this->_is_negative = true;
	else this->_is_negative = false;
	this->_digits.push_back(std::abs(i) % bigint::BASE);
	i /= bigint::BASE;
	if (i != 0) this->_digits.push_back(std::abs(i));
}

bigint::bigint(unsigned int i)
{
	this->_digits.push_back(i % bigint::BASE);
	i /= bigint::BASE;
	if (i != 0) this->_digits.push_back(i);
}

bigint::bigint(signed long l)
{
	if (l < 0) this->_is_negative = true;
	else this->_is_negative = false;
	this->_digits.push_back(std::abs(l) % bigint::BASE);
	l /= bigint::BASE;
	if (l != 0) this->_digits.push_back(std::abs(l));
}

bigint::bigint(unsigned long l)
{
	this->_digits.push_back(l % bigint::BASE);
	l /= bigint::BASE;
	if (l != 0) this->_digits.push_back(l);
}

bigint::bigint(signed long long l)
{
	if (l < 0) { this->_is_negative = true; l = -l; }
	else this->_is_negative = false;
	do {
		this->_digits.push_back(l % bigint::BASE);
		l /= bigint::BASE;
	} while (l != 0);
}

bigint::bigint(unsigned long long l)
{
	this->_is_negative = false;
	do {
		this->_digits.push_back(l % bigint::BASE);
		l /= bigint::BASE;
	} while (l != 0);
}

const bigint bigint::operator ++(int) // постфіксний інкремент
{
	*this += 1;
	return *this - 1;
}

const bigint bigint::operator --() // префіксний декремент
{
	return *this -= 1;
}

const bigint bigint::operator --(int) // постфіксний декремент
{
	*this -= 1;
	return *this + 1;
}

const bigint operator -(bigint left, const bigint& right) //віднімання двух чисел 
{

	if (right._is_negative) return left + (-right);
	else if (left._is_negative) return -(-left + right);
	else if (left < right) return -(right - left);
	int carry = 0;
	for (size_t i = 0; i < right._digits.size() || carry != 0; ++i)
	{
		left._digits[i] -= carry + (i < right._digits.size() ? right._digits[i] : 0);
		carry = left._digits[i] < 0;
		if (carry != 0) left._digits[i] += bigint::BASE;
	}

	left._remove_leading_zeros();
	return left;
}

bigint& bigint::operator -=(const bigint& value)
{
	return *this = (*this - value);
}

const bigint operator *(const bigint& left, const bigint& right) { // множення двух чисел
	bigint result;
	result._digits.resize(left._digits.size() + right._digits.size());
	for (size_t i = 0; i < left._digits.size(); ++i) {
		int carry = 0;
		for (size_t j = 0; j < right._digits.size() || carry != 0; ++j) {
			long long cur = result._digits[i + j] +
				left._digits[i] * 1LL * (j < right._digits.size() ? right._digits[j] : 0) + carry;
			result._digits[i + j] = static_cast<int>(cur % bigint::BASE);
			carry = static_cast<int>(cur / bigint::BASE);
		}
	}
	// не забудем про знак
	result._is_negative = left._is_negative != right._is_negative;
	result._remove_leading_zeros();
	return result;
}

bigint& bigint::operator *=(const bigint& value)
{
	return *this = (*this * value);
}

void bigint::_shift_right() // здвигає всі розряди на 1 вправо (домножує на BASE)
{
	if (this->_digits.size() == 0)
	{
		this->_digits.push_back(0);
		return;
	}
	this->_digits.push_back(this->_digits[this->_digits.size() - 1]);
	for (size_t i = this->_digits.size() - 2; i > 0; --i) this->_digits[i] = this->_digits[i - 1];
	this->_digits[0] = 0;
}

const bigint operator /(const bigint& left, const bigint& right) // ділення двух чисел
{
	if (right == 0) { cerr << "error" << endl; exit(0); }
	bigint b = right;
	b._is_negative = false;
	bigint result, current;
	result._digits.resize(left._digits.size());
	for (long long i = static_cast<long long>(left._digits.size()) - 1; i >= 0; --i) {
		current._shift_right();
		current._digits[0] = left._digits[i];
		current._remove_leading_zeros();
		int x = 0, l = 0, r = bigint::BASE;
		while (l <= r) {
			int m = (l + r) / 2;
			bigint t = b * m;
			if (t <= current) {
				x = m;
				l = m + 1;
			}
			else r = m - 1;
		}

		result._digits[i] = x;
		current = current - b * x;
	}

	result._is_negative = left._is_negative != right._is_negative;
	result._remove_leading_zeros();
	return result;
}

bigint& bigint::operator /=(const bigint& value)
{
	return *this = (*this / value);
}

const bigint operator %(const bigint& left, const bigint& right) // повертає решту від ділення
{
	bigint result = left - (left / right) * right;
	if (result._is_negative) result += right;
	return result;
}

bigint& bigint::operator %=(const bigint& value) // присвоює решту від ділення
{
	return *this = (*this % value);
}

bool bigint::odd() const // непарне
{
	if (this->_digits.size() == 0) return false;
	return this->_digits[0] & 1;
}

bool bigint::even() const // парне
{
	return !this->odd();
}

const bigint bigint::pow(bigint n) const // степінь
{
	bigint a(*this), result(1);
	while (n != 0)
	{
		if (n.odd()) result *= a;
		a *= a;
		n /= 2;
	}
	return result;
}

bigint bigint::sqrt() // корінь
{
	bigint l, r = *this;
	bigint res;
	while (l <= r)
	{
		bigint m = (l + r) / 2;
		if (m * m <= *this)
		{
			res = m;
			l = m + 1;
		}
		else
			r = m - 1;
	}
	return res;

}

bigint bigint::gcd(bigint left, bigint right) // нсд
{
	bigint nod = string("1");

	if (left > right) {
		bigint tmp = left;
		left = right;
		right = tmp;
	}

	while (left > long long(1) && right > long long(1)) {
		for (long long i = 2; i <= left; i++)
		{
			if (left % i == long long(0) && right % i == long long(0))
			{
				nod *= i;
				left /= i;
				right /= i;
				break;
			}
			if (left % i == long long(0))
			{
				left /= i;
				break;
			}
			if (right % i == long long(0))
			{
				right /= i;
				break;
			}
		}
	}
	return nod;
}

bigint bigint::abs(bigint x) // модуль
{
	if (x < 0) return -x;
	return x;
}

bigint fun_Elera(bigint mod)
{
	bigint result = mod;
	for (long long i = 2; i * i <= mod; ++i)
		if (mod % i == long long(0))
		{
			while (mod % i == long long(0))
				mod /= i;
			result = result - result / i;
		}
	if (mod > long long(1))
		result = result - result / mod;
	return result;
}

//https://www.mathemania.com/lesson/system-linear-congruences/
bigint bigint::congruencesOfTheFirstDegreeOne(bool* flagPtr, bigint a, bigint b, bigint mod, bigint* saveMod = 0)
{
	bigint result;
	bigint d = d.gcd(a, mod);
	if (d == 1)
	{
		result = b * a.pow(fun_Elera(mod) - 1) % mod;
	}
	else if (d > 1 && b % d == 0)  //d - рішень 
	{
		a = a / d;
		b = b / d;
		mod = mod / d;
		d = a % mod;
		if (d == 1)
		{
			if (saveMod != NULL) *saveMod = mod; // для підстановки
			result = b * a.pow(fun_Elera(mod) - 1) % mod;
		}
		else if (d > 1)
		{
			result = b * a.pow(fun_Elera(mod) - 1) % mod;
		}
		else
		{
			*flagPtr = 1;
			return bigint("-1");
		}

	}
	else
	{
		*flagPtr = 1;
		return bigint("-1");
	}
	return result;
}

bigint bigint::chineseRemainderTheorem(bigint* aB, bigint* bB, bigint* modB, size_t n, string* saveMe = 0) // тільки для взаємно простих чисел китайска теорема про решту
{
	bigint* x_b = new bigint[n];

	bigint* bBsave = new bigint[n];
	bigint* aBsave = new bigint[n];
	bigint* modBsave = new bigint[n];

	for (size_t i = 0; i < n; ++i)
	{
		bBsave[i] = bB[i];
		aBsave[i] = aB[i];
		modBsave[i] = modB[i];
	}

	bool flagPtr = 0;
	for (size_t i = 0; i < n; ++i)
	{
		x_b[i].congruencesOfTheFirstDegreeOne(&flagPtr, aB[i], bB[i], modB[i]);
		if (flagPtr)
		{
			*saveMe = "Система не має розв'язку. Тому що одне із рівнянь не вирішується!\n";
			return bigint("-1");
		}

	}

	bigint m = string("1");
	bigint* M = new bigint[n];
	bigint* M_ = new bigint[n];
	bigint result;

	for (size_t i = 0; i < n; ++i)
	{
		m = m * modB[i];
	}

	for (size_t i = 0; i < n; ++i)
	{

		M[i] = m / modB[i];
		M_[i] = M_[i].congruencesOfTheFirstDegreeOne(&flagPtr, M[i], 1, modB[i]);
		if (flagPtr)
		{
			*saveMe = "Система не має розв'язку. Тому що одне із рівнянь не вирішується!\n";
			return bigint("-1");
		}
		result = result + M[i] * M_[i] * bB[i];
	}
	result = result % m;

	*saveMe = "Система порівнянь\n\t || || ||\n\t \\/ \\/ \\/\n";
	for (size_t i = 0; i < n; ++i)
	{
		*saveMe += aBsave[i].toStr() + "*x = " + bBsave[i].toStr() + " + " + "(mod " + modBsave[i].toStr() + ')' + '\n';
	}
	*saveMe += "КИТАЙ\nМає розв'язок => x = " + result.toStr() + " + " + "(mod " + m.toStr() + ')' + "\n";

	return result;
}

bigint bigint::findMinVal(bigint* arr, size_t n, size_t tempI[])
{
	bigint min;
	min = arr[tempI[0]];//массива в переменные

	for (int r = 1; r < sizeof(tempI) / sizeof(size_t); r++)
	{
		if (min > arr[tempI[++r]]) min = arr[r];
	}
	return min;
}

bigint bigint::successiveSubstitutionMethod(bigint* aB, bigint* bB, bigint* modB, size_t n, string* saveMe = 0) // метод последовательной подстановки  
{
	bigint* x = new bigint[n];
	bigint g;
	bigint* M = new bigint[n];;

	bigint* bBsave = new bigint[n];
	bigint* aBsave = new bigint[n];
	bigint* modBsave = new bigint[n];

	string amm = "";

	for (size_t i = 0; i < n; ++i)
	{
		bBsave[i] = bB[i];
		aBsave[i] = aB[i];
		modBsave[i] = modB[i];
	}
	bool flagPtr = 0;
	bigint dd;

	for (size_t i = 0; i < n; ++i)
	{
		if (flagPtr)
		{
			*saveMe = "Система не має розв'язку!\n"
				"A imenno -> " + aBsave[i].toStr()
				+ "*x = " + bBsave[i].toStr()
				+ " + " + "(mod " + modBsave[i].toStr()
				+ ')' + '\n';
			return bigint("0");
		}
		x[i] = x[i].congruencesOfTheFirstDegreeOne(&flagPtr, aB[i], bB[i], modB[i]);
		bB[i] = x[i];
	}
	size_t k = 1;											//  ____________________________
															// |						    |
	for (size_t i = 0; i < n; ++i)							// | перевірка на вирішуваність |
	{														// | системи				    |
		for (; k < n; ++k)									// |____________________________|
		{
			dd = dd.gcd(modB[i], modB[k]);
			if ((bB[i] - bB[k]) % dd != long long(0))
			{
				*saveMe = "Sistema nie imejet reshenij.!\n";
				return bigint("-1");
			}
		}
		k = i + 2;
	}

	bigint d;
	bigint* tempG = new bigint[n];
	vector<size_t> tempI;

	size_t j = 1;
	for (size_t i = 0; i < n; ++i)
	{
		for (; j < n; ++j)
		{
			g = g.gcd(modB[i], modB[j]);
			if (g != long long(1))
			{
				tempG[i] = g;
				tempI.push_back((int)i);
			}
		}
		j = i + 2;
	}
	if (tempI.size() == 0)
	{
		tempI.push_back(1);
		d = long long(1);
	}
	else
	{
		size_t* arrI = new size_t[tempI.size()];

		for (size_t i = 0; i < tempI.size(); ++i)
		{
			arrI[i] = tempI[i];
		}
		d = d.findMinVal(tempG, n, arrI);
	}

	bigint AA_4et("1");
	bigint BB_4et("0");

	bigint AA_n4et("1");
	bigint MM_n4et;
	bigint BB_n4et("0");

	bigint tempA;
	bigint tempZ;
	bigint temp;
	bigint tempMod = modB[0];

	//первая итерация для чет
	AA_4et = modB[0];
	BB_4et = bB[1] - bB[0];
	tempA = tempA.congruencesOfTheFirstDegreeOne(&flagPtr, AA_4et, BB_4et, modB[0 + 1], &tempMod);
	//tempMod нужен для того, если в сравнении больше одного решения( d > 1 ), потому что там идёт
	//сокращение

	if (n > 2) tempMod = modB[2];

	//первая итерация для нечет
	AA_n4et = 1;
	BB_n4et = bB[0] + modB[0] * tempA;
	MM_n4et = modB[0] * modB[1];
	tempZ = tempA.congruencesOfTheFirstDegreeOne(&flagPtr, AA_n4et, BB_n4et, MM_n4et);
	temp = tempZ;

	//modB[0] = tempMod;
	for (size_t i = 1; i < n - 1; ++i) //    _____________________________
	{								   //	| решение методом подстановки |
		tempMod = modB[i + 1];		   //	|_____________________________|
		AA_4et = AA_4et * modB[i];
		BB_4et = bB[i + 1] - tempZ;

		tempA = tempA.congruencesOfTheFirstDegreeOne(&flagPtr, AA_4et, BB_4et, modB[i + 1], &tempMod);
		modB[i + 1] = tempMod;

		BB_n4et = temp + AA_4et * tempA;
		MM_n4et = MM_n4et * modB[i + 1];
		tempZ = tempA.congruencesOfTheFirstDegreeOne(&flagPtr, 1, BB_n4et, MM_n4et);
		temp = tempZ; // tempZ == b в конце

	}
	MM_n4et = MM_n4et / d; // 228
	*saveMe = "Sistema srawninij\n\t || || ||\n\t \\/ \\/ \\/\n";
	*saveMe = "Sistema srawninij\n\t || || ||\n\t \\/ \\/ \\/\n";
	for (size_t i = 0; i < n; ++i)
	{
		*saveMe += aBsave[i].toStr() + "*x = " + bBsave[i].toStr() +
			" + " + "(mod " + modBsave[i].toStr() + ')' + '\n';
	}

	*saveMe += "PODSTANOWKA\nImeet reshenije => x = " + tempZ.toStr() +
		" + " + "(mod " + MM_n4et.toStr() + ')' + "\n";
	return tempZ;
}

bigint bigint::solutionOfTheComparisonSystem(size_t n, bool showResult) // розв'язання систем порівнянь
{
	bigint result;

	bigint* x = new bigint[n];
	bigint g;
	bigint d;

	bigint* aB = new bigint[n];
	bigint* bB = new bigint[n];
	bigint* modB = new bigint[n];


	string* a = new string[n];
	string* b = new string[n];
	string* mod = new string[n];

	for (size_t i = 0; i < n; ++i)
	{
		cout << "Enter a >>> ";
		cin >> a[i];
		aB[i] = a[i];
		cout << "Enter b >>> ";
		cin >> b[i];
		bB[i] = b[i];
		cout << "Enter mod >>> ";
		cin >> mod[i];
		modB[i] = mod[i];
		cout << endl;
	}


	bool flagPtr = 0;


	if (n == 1)
	{
		result = result.congruencesOfTheFirstDegreeOne(&flagPtr, aB[0], bB[0], modB[0]);
		if (flagPtr)
		{
			cout << " nie imejet reshenij. Potomu 4to  nie reshajemo!\nKonkretno - >";
			cout << aB[0].toStr() + "*x = " + bB[0].toStr() + " + " + "(mod " + modB[0].toStr() + ')' << endl;
			result = long long(-1);
			return bigint(-1);
		}
		else
		{
			cout << aB[0].toStr() + "*x = " + bB[0].toStr() + " + " + "(mod " + modB[0].toStr() + ')' << endl;
			cout << "ELERA -> x = " + result.toStr() + " + " + "(mod " + modB[0].toStr() + ')' << endl;
			return result;
		}

	}
	else if (n > 1)
	{
		bigint* tempG = new bigint[n];
		vector<int> tempI;

		size_t j = 1;
		for (size_t i = 0; i < n; ++i)
		{
			for (; j < n; ++j)
			{
				g = g.gcd(modB[i], modB[j]);
				if (g != long long(1))
				{
					tempG[i] = g;
					tempI.push_back(i);
				}
			}
			j = i + 2;
		}

		if (tempI.size() > 0)
		{
			size_t* arrI = new size_t[tempI.size()];
			for (size_t i = 0; i < tempI.size(); ++i)
			{
				arrI[i] = tempI[i];
			}
			d = d.findMinVal(tempG, n, arrI);
		}
		else
			d = long long(1);

		if (d == 1)
		{
			string save;
			result = result.chineseRemainderTheorem(aB, bB, modB, n, &save);
			if (showResult) cout << save << endl;
			return result;
		}
		else
		{
			string save;
			auto start = chrono::high_resolution_clock::now();
			result = result.successiveSubstitutionMethod(aB, bB, modB, n, &save);
			auto end = chrono::high_resolution_clock::now();
			chrono::duration<float> duration = end - start;
			cout << "Duration = " << duration.count() << endl;
			if (showResult) cout << save << endl;
			return result;
		}
	}
	else
	{
		return bigint("-1");
	}
}

