package au.edu.uq.imb.memesuite.db;

import java.io.File;
import java.util.ArrayList;
import java.util.EnumSet;
import java.sql.*;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

import au.edu.uq.imb.memesuite.data.AlphStd;
import org.sqlite.*;

import static au.edu.uq.imb.memesuite.db.SQL.*;

public class DBList {
  private File db;
  private boolean deleteOnCleanup;
  protected SQLiteDataSource ds;
  private static Logger logger = Logger.getLogger("au.edu.uq.imb.memesuite.web");

  public DBList(File db, boolean deleteOnCleanup) throws ClassNotFoundException {
    Class.forName("org.sqlite.JDBC");
    SQLiteConfig config = new SQLiteConfig();
    config.enforceForeignKeys(true);
    config.setReadOnly(true);
    this.ds = new SQLiteDataSource(config);
    this.ds.setUrl("jdbc:sqlite:" + db);
    this.db = db;
    this.deleteOnCleanup = deleteOnCleanup;
  }

  public DBList(SQLiteDataSource ds) {
    this.ds = ds;
    this.db = null;
    this.deleteOnCleanup = false;
  }

  /**
   * Deletes the database file if the deleteOnCleanup flag is set.
   * Obviously you should only call this if you no-longer need the database!
   */
  public void cleanup() {
    if (deleteOnCleanup) {
      logger.log(Level.INFO, "Deleting unused database file \"" + this.db + "\".");
      if (!db.delete()) {
        logger.log(Level.WARNING, "Failed to delete database file \"" + this.db + "\".");
      }
    }
  }

  /**
   * Return the last modified date of the database.
   * @return the last modified date of the database.
   */
  public long lastModified() throws UnsupportedOperationException {
    if (db == null) throw new UnsupportedOperationException("DBList was initialised without a reference to the file so a last modified date can not be accessed.");
    return db.lastModified();
  }

  /**
   * Count the categories
   * @return the number of categories.
   * @throws SQLException if something goes wrong.
   */
  public int countCategories() throws SQLException {
    int count = 0;
    Connection connection = null;
    PreparedStatement ps = null;
    ResultSet rs = null;
    try {
      connection = ds.getConnection();
      ps = connection.prepareStatement(SQL.COUNT_CATEGORIES);
      rs = ps.executeQuery();
      if (!rs.next()) throw new SQLException("COUNT(*) should always return a row");
      count = rs.getInt(1);
      rs.close(); rs = null;
      ps.close(); ps = null;
      connection.close(); connection = null;
    } finally {
      if (rs != null) {
        try {
          rs.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (ps != null) {
        try {
          ps.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (connection != null) {
        try {
          connection.close();
        } catch (SQLException e) { /* ignore */ }
      }
    }
    return count;
  }
  protected PreparedStatement prepareCategoriesQuery(Connection conn,
      boolean shortOnly, Set<AlphStd> allowedAlphabets) throws SQLException {
    return conn.prepareStatement(SQL.SELECT_UNSORTED_CATEGORIES);
  }

  public List<Category> getCategories(boolean shortOnly, Set<AlphStd> allowedAlphabets) throws SQLException {
    CategoryView view = new CategoryView(shortOnly, allowedAlphabets);
    return view.items;
  }

  public int countListings() throws SQLException {
    int count = 0;
    Connection connection = null;
    PreparedStatement ps = null;
    ResultSet rs = null;
    try {
      connection = ds.getConnection();
      ps = connection.prepareStatement(SQL.COUNT_ALL_LISTINGS);
      rs = ps.executeQuery();
      if (!rs.next()) throw new SQLException("COUNT(*) should always return a row");
      count = rs.getInt(1);
      rs.close(); rs = null;
      ps.close(); ps = null;
      connection.close(); connection = null;
    } finally {
      if (rs != null) {
        try {
          rs.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (ps != null) {
        try {
          ps.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (connection != null) {
        try {
          connection.close();
        } catch (SQLException e) { /* ignore */ }
      }
    }
    return count;
  }

  public int countListings(long categoryId) throws SQLException {
    int count = 0;
    Connection connection = null;
    PreparedStatement ps = null;
    ResultSet rs = null;
    try {
      connection = ds.getConnection();
      ps = connection.prepareStatement(SQL.COUNT_CATEGORY_LISTINGS);
      ps.setLong(1, categoryId);
      rs = ps.executeQuery();
      if (!rs.next()) throw new SQLException("COUNT(*) should always return a row");
      count = rs.getInt(1);
      rs.close(); rs = null;
      ps.close(); ps = null;
      connection.close(); connection = null;
    } finally {
      if (rs != null) {
        try {
          rs.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (ps != null) {
        try {
          ps.close();
        } catch (SQLException e) { /* ignore */ }
      }
      if (connection != null) {
        try {
          connection.close();
        } catch (SQLException e) { /* ignore */ }
      }
    }
    return count;
  }

  protected PreparedStatement prepareListingsQuery(Connection conn, long categoryId,
      boolean shortOnly, EnumSet<AlphStd> allowedAlphabets) throws SQLException {
    PreparedStatement ps = conn.prepareStatement(SELECT_LISTINGS_OF_CATEGORY);
    ps.setLong(1, categoryId);
    return ps;
  }

  public List<Listing> getListings(long categoryId, boolean shortOnly,
      EnumSet<AlphStd> allowedAlphabets) throws SQLException {
    if (allowedAlphabets == null) throw new NullPointerException(
        "Parameter allowedAlphabets must not be null.");
    if (allowedAlphabets.isEmpty()) throw new IllegalArgumentException(
        "At least one allowed alphabet must be supplied.");
    ListingView view = new ListingView(categoryId, shortOnly, allowedAlphabets);
    return view.items;
  }

  public List<Listing> getListings(long categoryId) throws SQLException {
    ListingView view = new ListingView(categoryId, false, EnumSet.allOf(AlphStd.class));
    return view.items;
  }

  /**
   * Class for providing a list of objects retrieved from a
   * database
   */
  protected abstract class AbstractView<E> {
    protected List<E> items;

    /**
     * Construct a view.
     * Subclasses should call connect after doing their initialisation.
     */
    public AbstractView() {
      items = new ArrayList<E>();
    }

    /**
     * Run a query to create a result set for iterating over.
     */
    protected abstract PreparedStatement prepareStatement(Connection conn) throws SQLException;

    /**
     * Read values from a result set to construct a storage object.
     */
    protected abstract E constructData(ResultSet rset) throws SQLException;

    /**
     * Connect to the database and query the database for the list to iterate over.
     */
    protected void connect() throws SQLException {
      Connection conn = null;
      PreparedStatement stmt = null;
      ResultSet rset = null;
      try {
        conn = DBList.this.ds.getConnection();
        stmt = prepareStatement(conn);
        rset = stmt.executeQuery();
        while (rset.next()) {
          E item = constructData(rset);
          items.add(item);
        }
        rset.close();
        rset = null;
        stmt.close();
        stmt = null;
      } finally {
        if (rset != null) {
          try {
            rset.close();
          } catch (SQLException e) { /* ignore */ }
        }
        if (stmt != null) {
          try {
            stmt.close();
          } catch (SQLException e) { /* ignore */ }
        }
        if (conn != null) {
          try {
            conn.close();
          } catch (SQLException e) { /* ignore */ }
        }
      }
    }
  }

  public class CategoryView extends AbstractView<Category> {
    protected boolean shortOnly;
    protected Set<AlphStd> allowedAlphabets;

    public CategoryView(boolean shortOnly, Set<AlphStd> allowedAlphabets) throws SQLException {
      this.shortOnly = shortOnly;
      this.allowedAlphabets = allowedAlphabets;
      connect();
    }

    protected PreparedStatement prepareStatement(Connection conn) throws SQLException {
      return prepareCategoriesQuery(conn, shortOnly, allowedAlphabets);
    }

    protected Category constructData(ResultSet rset) throws SQLException {
      long id = rset.getLong(1);
      String name = rset.getString(2);
      long count = rset.getLong(3);
      boolean priors = rset.getBoolean(4);
      return new Category(id, name, count, priors);
    }
  }

  public class ListingView extends AbstractView<Listing> {
    protected long categoryId;
    protected boolean shortOnly;
    protected EnumSet<AlphStd> allowedAlphabets;

    public ListingView(long categoryId, boolean shortOnly, EnumSet<AlphStd> allowedAlphabets) throws SQLException {
      this.categoryId = categoryId;
      this.shortOnly = shortOnly;
      this.allowedAlphabets = allowedAlphabets;
      connect();
    }

    protected PreparedStatement prepareStatement(Connection conn) throws SQLException {
      return prepareListingsQuery(conn, categoryId, shortOnly, allowedAlphabets);
    }

    protected Listing constructData(ResultSet rset) throws SQLException {
      long id = rset.getLong(1);
      String name = rset.getString(2);
      String description = rset.getString(3);
      EnumSet<AlphStd> alphabets = SQL.intToEnums(AlphStd.class, rset.getInt(4));
      boolean priors = rset.getBoolean(5);
      return new Listing(id, name, description, alphabets, priors);
    }
  }


}
