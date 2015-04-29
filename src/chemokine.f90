! Chemokine data

module chemokine

use global

implicit none

type receptor_type
	character(12) :: name
	logical :: used
	integer :: chemokine
	integer :: sign							! 1 = attraction, -1 = repulsion (?)
	real(REAL_KIND) :: strength				! relative strength of chemotactic influence
!	real(REAL_KIND) :: level(5)				! receptor levels for different cell stages
	real(REAL_KIND) :: saturation_threshold	! chemokine signal(gradient) level producing desensitization
	real(REAL_KIND) :: refractory_time		! duration of desensitization
end type

type chemokine_type
	character(12) :: name
	logical :: used
	logical :: use_secretion
	real(REAL_KIND) :: bdry_rate
	real(REAL_KIND) :: bdry_conc
	real(REAL_KIND) :: diff_coef
	real(REAL_KIND) :: halflife
	real(REAL_KIND) :: decay_rate
	real(REAL_KIND) :: radius
	real(REAL_KIND), allocatable :: coef(:,:)
	real(REAL_KIND), allocatable :: conc(:,:,:)
	real(REAL_KIND), allocatable :: grad(:,:,:,:)
end type

type ODEdiff_type
	integer :: ichemo
	integer :: nvars
	integer, allocatable :: ivar(:,:,:)
	integer, allocatable :: varsite(:,:)
	integer, allocatable :: icoef(:,:)
end type

type(chemokine_type), target :: chemo(MAX_CHEMO)
type(receptor_type), target :: receptor(MAX_RECEPTOR)
type(ODEdiff_type) :: ODEdiff

end module