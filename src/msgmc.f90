module msgmc
  interface
    subroutine msgmci (i, a)
      integer i, a
    end subroutine msgmci
  end interface

  interface
    subroutine msgmca
    end subroutine msgmca
  end interface

  interface
    subroutine msgmce (a)
      integer a
    end subroutine msgmce
  end interface

  interface
    subroutine msgmcl
    end subroutine msgmcl
  end interface
end module msgmc
