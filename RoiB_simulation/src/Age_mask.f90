SUBROUTINE Age_mask( Model,Solver,dt,TransientSimulation )
    USE DefUtils
    IMPLICIT NONE
    !------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    TYPE(Model_t) :: Model
    REAL(KIND=dp) :: dt
    LOGICAL :: TransientSimulation
    !------------------------------------------------------------------------------
    ! Local variables
    !------------------------------------------------------------------------------
    Integer                       :: DIM,N,e,i,j,k,kf,km,ka,kd,kd2,kk,NMax,NM,stat,NMB,NB
    Logical                       :: GotIt,FirstTime=.TRUE.,IsThere,ALE
    REAL(KIND=dp), POINTER        :: VariableValues(:)
    TYPE(Element_t),POINTER       :: Element, UElement, ElementB
    TYPE(Nodes_t)                 :: ElementNodes
    INTEGER,  POINTER             :: NodeIndexes(:), NodeIndexesB(:)
    REAL(KIND=dp),ALLOCATABLE,SAVE:: Basis(:), dBasisdx(:,:),ddBasisddx(:,:,:)
    TYPE(Matrix_t), POINTER       :: Systemmatrix
    CHARACTER(LEN=MAX_NAME_LEN)   :: FlowSolName
    TYPE(Variable_t), POINTER     :: FlowSol,MeshVelSol, PointerToVariable
    INTEGER, POINTER              :: FlowPerm(:),MeshVelPerm(:)
    REAL(KIND=dp), POINTER        :: Flow(:),MeshVel(:)
    INTEGER                       :: FlowDOFs,MeshVelDOFs
    INTEGER                       :: NNMAX,NNMAXB
    REAL(KIND=dp)                 :: xa(3),xm(3),xd(3),alpha(3),um(3)
    REAL(KIND=dp), DIMENSION(3)   :: LocalCoordinates
    REAL(KIND=dp)                 :: eps=1.0e-6, eps2=1.0e-6
    INTEGER                       :: MaxDOFs
    REAL(KIND=dp), ALLOCATABLE    :: Vector(:,:)
    REAL(KIND=dp)                 :: Agea,RcvAgea,TmpAge,BMBa,RcvBMBa,TmpBMB
    REAL(KIND=dp)                 :: s1(3),s2(3),PoR(3),d2,d2min
    INTEGER                       :: emin
    REAL(KIND=dp)                 :: at,st,totat,totst,Norm,PrevNorm,RelativeChange
    LOGICAL, ALLOCATABLE          :: Found(:),isBC(:),isAgeIn(:),isAtBoundary(:)
    LOGICAL, ALLOCATABLE          :: SkipThisNodeInThisPartition(:)
    LOGICAL                       :: Owner=.FALSE., AtInterface=.FALSE., UnFoundFatal=.TRUE.
    REAL(KIND=dp), ALLOCATABLE    :: Cond(:)
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER                       :: nn,precv
    INTEGER, POINTER              :: nlist(:)
    INTEGER                       :: ierr,gk
    INTEGER                       :: request(ParEnv % PEs)
    TYPE buffer_t
       INTEGER                    :: n
       INTEGER, ALLOCATABLE       :: gbuff(:)
       REAL(KIND=dp), ALLOCATABLE :: vbuff(:)
    END TYPE buffer_t
    TYPE(buffer_t)                :: RequestSend(ParEnv % PEs),RequestRecv(ParEnv % PEs)
    TYPE(buffer_t)                :: ReplySend(ParEnv % PEs),ReplyRecv(ParEnv % PEs)
!!!!!!!!!!!!!!!!!!! [VV] !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TYPE(ValueList_t), POINTER    :: SolverParams, BodyForce, BC, BMB2
    TYPE(Variable_t), POINTER     :: Age,GM
    REAL(KIND=dp), POINTER        :: Age_val(:),GM_val(:)
    INTEGER, POINTER              :: Age_Perm(:),GM_Perm(:)
    REAL(KIND=dp)                 :: rhoi, rhoo, yearinsec
    REAL(KIND=dp),SAVE            :: Lambda
    real(kind=dp)                 :: Coord(3),UVW(3)
    INTEGER,SAVE                  :: VDOFs
    REAL(KIND=dp), ALLOCATABLE    :: BMB1(:), SMB1(:)
    INTEGER                       :: bf_id, istat, NodeI_m
    REAL(KIND=dp)                 :: u,v,w,SqrtElementMetric,s
    TYPE(Nodes_t)                 :: Nodes, NodesB
    TYPE(GaussIntegrationPoints_t):: IP
    REAL(KIND=dp), ALLOCATABLE    :: MASS(:,:), STIFF(:,:), LOADB(:),LOADS(:), FORCE(:)
    Real(KIND=dp)                 :: xx(3),yy(3),yymin(3),dist2,dist2min
    REAL(KIND=dp)                 :: Vza
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(Message,'(a)') 'Start Solver'
    CALL Info('Vz_SSA_BMB_SMB', Message, Level=4)
    !------------------------------------------------------------------------------
    ! Get Constants
    !------------------------------------------------------------------------------
    DIM           = CoordinateSystemDimension()
    NMAX          = Model  % MaxElementNodes
    NM            = Solver % Mesh % NumberOfNodes
    NMB           = GetNofBoundaryElements()
    ! ALLOCATE( FORCE(NMAX), LOADB(NMAX), LOADS(NMAX), MASS(NMAX,NMAX), STIFF(NMAX,NMAX), STAT=istat )
    IF(DIM        == 2)  THEN
       NNMAX      = 10!4
    ELSE
       NNMAX      = 20!8
    END IF
    MaxDOFs=DIM
    !------------------------------------------------------------------------------
    !    Get variables for the solution
    !------------------------------------------------------------------------------

    SolverParams => GetSolverParams()
    ! BodyForce    => GetBodyForce()
    BC           => GetBC()

    yearinsec    =  365.25*24*60*60
    rhoi         =  910.0/(1.0e6*yearinsec*yearinsec)
    rhoo         =  1025.0/(1.0e6*yearinsec*yearinsec)
    ! rhoi         =  GetConstReal( Model % Constants,'rhoi') ! Model % Constants
    ! rhoo         =  GetConstReal( Model % Constants,'rhoo')

    Age          => VariableGet( Model % Mesh % Variables, 'Age',UnFoundFatal=UnFoundFatal)
    Age_val      => Age  % Values
    Age_Perm     => Age  % Perm!
    !
    GM           => VariableGet( Model % Mesh % Variables, 'Age_mask',UnFoundFatal=UnFoundFatal)
    GM_val       => GM % Values
    GM_Perm      => GM % Perm!

    ALLOCATE( ElementNodes % x( NMAX ),         &
         ElementNodes % y( NMAX ),              &
         ElementNodes % z( NMAX ),              &
         Vector(MaxDOFS,NMAX),                  &
         BMB1(NMAX),                            &
         SMB1(NMAX),                            &
         ddBasisddx(n,3,3),                     &
         STAT=stat  )
     ALLOCATE(isBC(NM),isAtBoundary(NM),        &
         STAT=stat  )
     ALLOCATE(isAgeIn(NM),                      &
         STAT=stat  )
DO e = 1, Solver % NumberOfActiveElements ! main loop through the active elements
           Element             => GetActiveElement(e)
           NodeIndexes         => Element % NodeIndexes
           n                   =  GetElementNOFNodes(Element)
 Do kk = 1,n  !! loop through the nodes of the element (e)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        IF (GM_val(GM_Perm(NodeIndexes(kk)))  < 1.0 ) THEN
           Age_val(Age_Perm(NodeIndexes(kk))) = 500.0
        END IF
 End DO !! End of loop through the nodes n of the element (e)

END DO !! End of loop through all active elements - main loop


CONTAINS
    FUNCTION DistanceSquared(p1,p2) RESULT(d2)
      !Square of the distance between 2 points P1 and P2
      IMPLICIT NONE
      DOUBLE PRECISION :: p1(:),p2(:),d2

      d2=(p2(1)-p1(1))*(p2(1)-p1(1))+(p2(2)-p1(2))*(p2(2)-p1(2))+(p2(3)-p1(3))*(p2(3)-p1(3))

    END FUNCTION DistanceSquared


  End SUBROUTINE Age_mask
