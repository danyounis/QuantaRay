module prec
integer, parameter, public :: sngl_t = selected_real_kind(6,37)
integer, parameter, public :: dble_t = selected_real_kind(15,307)
#ifdef SINGLE
integer, parameter, public :: num = sngl_t
#elif DOUBLE
integer, parameter, public :: num = dble_t
#endif
end module prec
