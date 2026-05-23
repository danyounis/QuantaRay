module prec
integer, parameter, public :: sngl_t = selected_real_kind(6,37)
integer, parameter, public :: dble_t = selected_real_kind(15,307)
integer, parameter, public :: quad_t = selected_real_kind(32)
#ifdef SINGLE
integer, parameter, public :: num = sngl_t
#elif DOUBLE
integer, parameter, public :: num = dble_t
#elif QUADRUPLE
integer, parameter, public :: num = quad_t
#endif
end module prec
