//
// Created by alexander on 29.01.21.
//

#ifndef PAPILO_EARLYPRESOLVETERMINATIONEXCEPTION_HPP
#define PAPILO_EARLYPRESOLVETERMINATIONEXCEPTION_HPP

#include "papilo/core/PresolveMethod.hpp"

namespace papilo
{
class EarlyPresolveTerminationException : public std::exception
{

 private:
   const PresolveStatus presolve_status;

 public:
   explicit EarlyPresolveTerminationException( const PresolveStatus status ): presolve_status(status) {
   }

   PresolveStatus get_presolve_status() const { return presolve_status; }
};
} // namespace papilo


#endif // PAPILO_EARLYPRESOLVETERMINATIONEXCEPTION_HPP
