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

  interface
    subroutine msgmci2 (i, a, b)
      integer i, a, b
    end subroutine msgmci2
  end interface

  interface
    subroutine msgmca2
    end subroutine msgmca2
  end interface

  interface
    subroutine msgmce2 (a, b)
      integer a, b
    end subroutine msgmce2
  end interface

  interface
    subroutine msgmcl2
    end subroutine msgmcl2
  end interface
end module msgmc
