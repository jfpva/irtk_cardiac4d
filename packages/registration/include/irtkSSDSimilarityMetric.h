/*=========================================================================

  Library   : Image Registration Toolkit (IRTK)
  Module    : $Id$
  Copyright : Imperial College, Department of Computing
              Visual Information Processing (VIP), 2008 onwards
  Date      : $Date$
  Version   : $Revision$
  Changes   : $Author$

=========================================================================*/

#ifndef _IRTKSSDSIMILARITYMETRIC_H

#define _IRTKSSDSIMILARITYMETRIC_H

/**
 * Class for voxel similarity measures based on sums-of-squared differences.
 *
 */

class irtkSSDSimilarityMetric : public irtkSimilarityMetric
{

private:

  /// Number of samples
  double _n;

  /// Sums-of-squared differences
  double _ssd;

public:

  /// Constructor
  irtkSSDSimilarityMetric();

  /// Add sample
  virtual void Add(int, int);

  /// Remove sample
  virtual void Delete(int, int);

  /// Combine similarity metrics
  virtual void Combine(irtkSimilarityMetric *);

  /// Reset similarity metric
  virtual void Reset();

  /// Reset similarity metric
  virtual void Reset(irtkSimilarityMetric *);

  /// Evaluate similarity measure
  virtual double Evaluate();

  virtual int** GetBins();

  virtual void PutBins(int**);

  virtual void PutNSamp(int);

  virtual void PutSSDVal(int);
};

inline int** irtkSSDSimilarityMetric::GetBins()
{
	return NULL;
}

inline void irtkSSDSimilarityMetric::PutBins(int** _bins)
{
	_ssd = **_bins;
}

inline void irtkSSDSimilarityMetric::PutNSamp(int _nsamp)
{
	_n = _nsamp;
}

inline void irtkSSDSimilarityMetric::PutSSDVal(int ssd)
{
	_ssd = ssd;
}

inline irtkSSDSimilarityMetric::irtkSSDSimilarityMetric()
{
  _ssd = 0;
  _n = 0;
}

inline void irtkSSDSimilarityMetric::Add(int x, int y)
{
  _ssd += (x-y)*(x-y);
  _n++;
}

inline void irtkSSDSimilarityMetric::Delete(int x, int y)
{
  _ssd -= (x-y)*(x-y);
  _n--;
}

inline void irtkSSDSimilarityMetric::Combine(irtkSimilarityMetric *metric)
{
  irtkSSDSimilarityMetric *m = dynamic_cast<irtkSSDSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkSSDSimilarityMetric::Combine: Dynamic cast failed" << endl;
    exit(1);
  }

  _n   += m->_n;
  _ssd += m->_ssd;
}


inline void irtkSSDSimilarityMetric::Reset()
{
  _ssd = 0;
  _n = 0;
}

inline void irtkSSDSimilarityMetric::Reset(irtkSimilarityMetric *metric)
{
  irtkSSDSimilarityMetric *m = dynamic_cast<irtkSSDSimilarityMetric *>(metric);

  if (m == NULL) {
    cerr << "irtkSSDSimilarityMetric::Reset: Dynamic cast failed" << endl;
    exit(1);
  }

  _ssd = m->_ssd;
  _n = m->_n;
}

inline double irtkSSDSimilarityMetric::Evaluate()
{
  if (_n > 0) {
    return - _ssd / _n;
  } else {
    cerr << "irtkSSDSimilarityMetric::Evaluate: No samples";
    return 0;
  }
}

#endif
