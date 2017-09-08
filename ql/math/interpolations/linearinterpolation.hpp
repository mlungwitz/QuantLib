/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2000, 2001, 2002, 2003 RiskMap srl
 Copyright (C) 2003, 2004, 2008 StatPro Italia srl
 Copyright (C) 2017 Matthias Lungwitz

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file linearinterpolation.hpp
    \brief linear interpolation between discrete points
*/

#ifndef quantlib_linear_interpolation_hpp
#define quantlib_linear_interpolation_hpp

#include <ql/math/interpolation.hpp>
#include <vector>

namespace QuantLib {

    namespace detail {
        template<class I1, class I2> class LinearInterpolationImpl;
    }

    //! %Linear interpolation between discrete points
    /*! \ingroup interpolations */
    class LinearInterpolation : public Interpolation {
      public:
		  enum ExtrapolationType {
			  Flat,
			  Linear
		};

        /*! \pre the \f$ x \f$ values must be sorted. */
        template <class I1, class I2>
        LinearInterpolation(const I1& xBegin, const I1& xEnd,
                            const I2& yBegin,
							LinearInterpolation::ExtrapolationType leftCond = Linear,
							LinearInterpolation::ExtrapolationType rightCond = Linear) {
            impl_ = boost::shared_ptr<Interpolation::Impl>(new
                detail::LinearInterpolationImpl<I1,I2>(xBegin, xEnd,
                                                       yBegin, leftCond, rightCond));
            impl_->update();
        }
    };

    //! %Linear-interpolation factory and traits
    /*! \ingroup interpolations */
    class Linear {
      public:
        template <class I1, class I2>
        Interpolation interpolate(const I1& xBegin, const I1& xEnd,
                                  const I2& yBegin,
								  LinearInterpolation::ExtrapolationType leftCond = LinearInterpolation::Linear,
								  LinearInterpolation::ExtrapolationType rightCond = LinearInterpolation::Linear) const {
            return LinearInterpolation(xBegin, xEnd, yBegin, leftCond, rightCond);
        }
        static const bool global = false;
        static const Size requiredPoints = 2;
    };

    namespace detail {

        template <class I1, class I2>
        class LinearInterpolationImpl
            : public Interpolation::templateImpl<I1,I2> {
          public:
            LinearInterpolationImpl(const I1& xBegin, const I1& xEnd,
                                    const I2& yBegin, LinearInterpolation::ExtrapolationType leftCondition = Linear,
									LinearInterpolation::ExtrapolationType rightCondition = Linear)
            : Interpolation::templateImpl<I1,I2>(xBegin, xEnd, yBegin,
                                                 Linear::requiredPoints),
              primitiveConst_(xEnd-xBegin), s_(xEnd-xBegin),
		      leftType_(leftCondition), rightType_(rightCondition) {}
            void update() {
                primitiveConst_[0] = 0.0;
                for (Size i=1; i<Size(this->xEnd_-this->xBegin_); ++i) {
                    Real dx = this->xBegin_[i]-this->xBegin_[i-1];
                    s_[i-1] = (this->yBegin_[i]-this->yBegin_[i-1])/dx;
                    primitiveConst_[i] = primitiveConst_[i-1]
                        + dx*(this->yBegin_[i-1] +0.5*dx*s_[i-1]);
                }

				if (rightType_ == LinearInterpolation::Flat)
					primitiveConst_[Size(this->xEnd_ - this->xBegin_) - 1] = primitiveConst_[Size(this->xEnd_ - this->xBegin_) - 2];
            }
            Real value(Real x) const {
				Size i = this->locate(x);
				if ((leftType_ == LinearInterpolation::Flat) && x < xMin() )
					return this->yBegin_[0];
				else if ((rightType_ == LinearInterpolation::Flat) && x > xMax())
					return this->yBegin_[i];
				else
					return this->yBegin_[i] + (x-this->xBegin_[i])*s_[i];
            }
            Real primitive(Real x) const {
				Size i = this->locate(x);
				if ((leftType_ == LinearInterpolation::Flat) && x < xMin())
					return this->yBegin_[0] * (x - xMin());
				else if ((rightType_ == LinearInterpolation::Flat) && x > xMax())
				{
					Real dx1 = xMax() - this->xBegin_[i];
					Real dx2 = x - xMax();
					return primitiveConst_[i] +
						dx1*(this->yBegin_[i] + 0.5*dx1*s_[i]) +
						dx2*(this->yBegin_[i] + dx1*s_[i]);
				}		
				else
				{
					Real dx = x - this->xBegin_[i];
					return primitiveConst_[i] +
						dx*(this->yBegin_[i] + 0.5*dx*s_[i]);
				}

            }
            Real derivative(Real x) const {
				if ((leftType_ == LinearInterpolation::Flat) && x < xMin())
					return 0.0;
				else if ((rightType_ == LinearInterpolation::Flat) && x > xMax())
					return 0.0;
				else
				{
					Size i = this->locate(x);
					return s_[i];
				}
            }
            Real secondDerivative(Real) const {
                return 0.0;
            }
          private:
            std::vector<Real> primitiveConst_, s_;
			LinearInterpolation::ExtrapolationType leftType_, rightType_;
        };

    }

}

#endif
