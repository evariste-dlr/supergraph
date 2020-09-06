#ifndef __ECMAPPING_H__
#define __ECMAPPING_H__

/**
 * @file ECMapping.h
 *
 *  Defines the class ECMapping.
 *
 * @author Évariste Daller <<daller@unicaen.fr>>
 * @date 15/05/2018
 */


/**
 * @brief An Error-correcting mapping
 */
class ECMapping {

protected:
  
  int* _forward;  //< G1 -> G2
  int* _reverse;  //< G2 -> G1
  int _n;
  int _m;

  double _cost;

public:

  /**
   * @brief Construct a ECMapping from 2 existing arrays
   *
   *   Arrays aren't copied, but only referenced. This means that
   *   a particular attention should be given on shared memory issues.
   *   If you need the instance to let your data unharmed, use the const
   *   version instead.
   *
   * @param f  First set to second set mapping
   * @param r  Second set to first set mapping
   * @param n  Size of the first set
   * @param m  Size of the second set
   * @param cost  Cost of the mapping
   */
  ECMapping (int* f, int* r, int n, int m, double cost) :
    _forward(f),
    _reverse(r),
    _n(n), _m(m),
    _cost(cost)
  {}


  /**
   * @brief Construct a ECMapping from 2 existing arrays (const version)
   *
   *   Arrays are copied locally to the instance.
   *
   * @param f  First set to second set mapping
   * @param r  Second set to first set mapping
   * @param n  Size of the first set
   * @param m  Size of the second set
   * @param cost  Cost of the mapping
   */
  ECMapping (const int* f, const int* r, int n, int m, double cost) :
    _n(n), _m(m),
    _cost(cost)
  {
    _forward = new int[_n];
    _reverse = new int[_m];
    memcpy(_forward, f, (_n)*sizeof(int));
    memcpy(_reverse, r, (_m)*sizeof(int));
  }


  /**
   * @brief Copy constructor, performs memory copies of the arrays
   */
  ECMapping (const ECMapping & other) :
    _n(other._n),
    _m(other._m),
    _cost(other._cost)
  {
    _forward = new int[_n];
    _reverse = new int[_m];
    memcpy(_forward, other._forward, (_n)*sizeof(int));
    memcpy(_reverse, other._reverse, (_m)*sizeof(int));
  }


  ~ECMapping (){
    delete [] _forward;
    delete [] _reverse;
  }


  /**
   * @brief array subscripting : access the forward mapping
   * @see f(int)
   * @see r(int)
   */
  int& operator[] (size_t i){ return _forward[i]; }
  const int& operator[] (size_t i) const { return _forward[i]; }

  int* f () { return _forward; }
  int* r () { return _reverse; }

  int f (int i) const { return _forward[i]; }
  int r (int i) const { return _reverse[i]; }

  int n () const { return _n; } //!< Size of the first set
  int m () const { return _m; } //!< Size of the second set

  /**
   * @brief Cost of the mapping
   */
  double cost() const { return _cost; }


  /**
   * @brief  True iff operation attached to element `i<n` of the first set is a removal
   */
  bool is_rem(int i) const { return _forward[i] >= _m; }

  /**
   * @brief  True iff operation attached to element `j<m` of the second set is an insertion
   */
  bool is_ins(int j) const { return _reverse[j] >= _n; }



  /*****************************
   * Comparison functions
   *****************************/

  /**
   * @brief Lower or equal in term of cost
   */
  static bool leq( const ECMapping & m1, const ECMapping & m2 ){
    return m1.cost() <= m2.cost();
  }
  
  /**
   * @brief Lower than, in term of cost
   */
  static bool lt( const ECMapping & m1, const ECMapping & m2 ){
    return m1.cost() < m2.cost();
  }

  /**
   * @brief Greater or equal in term of cost
   */
  static bool geq( const ECMapping & m1, const ECMapping & m2 ){
    return m1.cost() >= m2.cost();
  }

  /**
   * @brief Greater than, in term of cost
   */
  static bool gt( const ECMapping & m1, const ECMapping & m2 ){
    return m1.cost() > m2.cost();
  }
  

  /**
   * @brief Lower or equal in term of cost (pointer version)
   */
  static bool leq_ptr( ECMapping* const & m1, ECMapping* const & m2 ){
    return m1->cost() <= m2->cost();
  }

  /**
   * @brief Lower than, in term of cost (pointer version)
   */
  static bool lt_ptr( ECMapping* const & m1, ECMapping* const & m2 ){
    return m1->cost() < m2->cost();
  }

  /**
   * @brief Greater or equal in term of cost (pointer version)
   */
  static bool geq_ptr( ECMapping* const & m1, ECMapping* const & m2 ){
    return m1->cost() >= m2->cost();
  }

  /**
   * @brief Greater than, in term of cost (pointer version)
   */
  static bool gt_ptr( ECMapping* const & m1, ECMapping* const & m2 ){
    return m1->cost() > m2->cost();
  }
};


#endif
