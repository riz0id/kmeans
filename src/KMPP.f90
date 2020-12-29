module kmpp

implicit none

contains

pure elemental function dist(x, y) result(d)
  real*8, intent(in) :: x
  real*8, intent(in) :: y
  real*8             :: d
  d = sqrt((x - y)**2)
end function dist

subroutine recluster(X, n, k, coids, distances)
  real*8,  dimension(n),    intent(in)    :: X
  integer,                  intent(in)    :: n
  integer,                  intent(in)    :: k
  real*8,  dimension(k),    intent(inout) :: coids
  real*8,  dimension(k, n), intent(inout) :: distances

  integer, dimension(k) :: lens
  integer :: i, j, min_ix
  real*8 :: min_dist, sum
  do i = 1, n
    min_dist = HUGE(min_dist)
    min_ix   = 1

    do j = 1, k
      if (distances(j, i) .LT. min_dist) then
        min_dist = distances(j, i)
        min_ix = j
      end if
    end do

    do j = 1, k
      if (j .NE. min_ix) then
        distances(j, i) = 0
      end if
    end do

    lens(min_ix) = lens(min_ix) + 1
  end do

  do i = 1, k
    sum = 0

    do j = 1, n
      sum = sum + distances(i, j)
    end do

    coids(k) = sum / lens(k)
  end do
end subroutine recluster

subroutine distances(X, n, k, coids, clusters)
  real*8,  dimension(n),    intent(in)    :: X
  integer,                  intent(in)    :: n
  integer,                  intent(in)    :: k
  real*8,  dimension(k),    intent(in)    :: coids
  real*8,  dimension(k, n), intent(inout) :: clusters

  integer :: i, j
  do i = 1, k
    do j = 1, n
      clusters(i, j) = dist(coids(i), X(j))
    end do
  end do
end subroutine distances

subroutine kmeans(X, n, k, coids, iters)
  real*8,  dimension(n), intent(in)    :: X
  integer,               intent(in)    :: n
  integer,               intent(in)    :: k
  real*8,  dimension(k), intent(inout) :: coids
  integer,               intent(in)    :: iters

  real*8, dimension(k, n) :: clusters

  integer :: iter
  do iter = 1, iters
    call distances(X, n, k, coids, clusters)
    call recluster(X, n, k, coids, clusters)
  end do
end subroutine kmeans

end module kmpp
